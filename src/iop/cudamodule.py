"""
CUDA Source Manager

This module provides a robust solution for managing CUDA source files in Python applications using CuPy.
It addresses several key challenges in CUDA development:

Key Features:
1. Preprocessor for Include Directives: Handles #include directives in CUDA source files, allowing 
   modular code organization similar to standard C/C++ development.

2. Line Number Mapping: Maintains a mapping between lines in the processed source and the original 
   files, enabling accurate error reporting.

3. Error Reporting Enhancement: Transforms cryptic CUDA compilation errors into readable messages 
   that point to the exact location in your original source files.

4. Circular Include Detection: Prevents infinite recursion by detecting circular dependencies in include files.

5. Source Caching: Improves performance by caching processed source files.

Usage:
The manager simplifies CUDA kernel development by:
1. Allowing code to be split across multiple files with #include directives
2. Providing clear error messages that reference your original source files
3. Showing available kernel functions when a function name is mistyped
4. Offering debugging tools like saving the processed source

This module bridges the gap between CuPy's raw CUDA capabilities and a more developer-friendly workflow,
making CUDA kernel development and debugging significantly easier in Python applications.
"""

import os   
import re
import cupy as cp
import sys

import iop.cudaKernelWrapper as ckw
from python.schism_consts import DEBUG

# ANSI color codes for terminal output
class Colors:
    RESET = "\033[0m"
    RED = "\033[91m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    BLUE = "\033[94m"
    MAGENTA = "\033[95m"
    CYAN = "\033[96m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"

class CudaSourceManager:
    def __init__(self, project_root):
        """Initialize the CUDA source manager with project root directory."""
        self.project_root = project_root
        self.source_cache = {}  # Cache for processed source files
        self.line_map = {}      # Maps combined file line numbers to original files and lines
        
    def preprocess_file(self, filepath, included_files=None, line_counter=None):
        """
        Preprocess a CUDA source file, handling #include directives.
        Returns the processed source and updates line mapping.
        
        Args:
            filepath: Path to the file to preprocess
            included_files: Set of files being processed (to detect circular includes)
            line_counter: Dictionary to track current line number in output
        """
        if included_files is None:
            included_files = set()
        
        if line_counter is None:
            line_counter = {"current": 1}
        
        # Check if file is already processed to avoid circular includes
        if filepath in included_files:
            return f"// Circular include detected: {filepath}\n"
        
        # Add to processed files
        included_files.add(filepath)
        
        # Check if file is in cache
        if filepath in self.source_cache:
            return self.source_cache[filepath]
        
        # Get directory of current file for relative includes
        file_dir = os.path.dirname(filepath)
        
        # Process the file
        processed_source = []
        
        with open(filepath, 'r', encoding='utf-8') as f:
            file_lines = f.readlines()
            
            for line_num, line in enumerate(file_lines, 1):
                # Check for include directive
                include_match = re.match(r'#include\s+"([^"]+)"', line)
                
                if include_match:
                    # Handle include directive
                    include_file = include_match.group(1)
                    include_path = os.path.normpath(os.path.join(file_dir, include_file))
                    
                    # Add a comment to mark the include
                    processed_source.append(f"// BEGIN INCLUDE: {include_file}")
                    self.line_map[line_counter["current"]] = (filepath, line_num, None)
                    line_counter["current"] += 1
                    
                    # Process the included file
                    included_source = self.preprocess_file(include_path, included_files, line_counter)
                    processed_source.append(included_source)
                    
                    # Add a comment to mark the end of include
                    processed_source.append(f"// END INCLUDE: {include_file}")
                    self.line_map[line_counter["current"]] = (filepath, line_num, None)
                    line_counter["current"] += 1
                else:
                    # Regular line, add to processed source
                    processed_source.append(line.rstrip('\n'))
                    
                    # Update line mapping - map the output line to the original file and line
                    self.line_map[line_counter["current"]] = (filepath, line_num, None)
                    line_counter["current"] += 1
        
        # Remove from processed files
        included_files.remove(filepath)
        
        # Join the processed source
        result = '\n'.join(processed_source)
        
        # Cache the result
        self.source_cache[filepath] = result
        
        return result
    
    
    def create_module(self, main_file, options=None):
        """
        Create a CuPy RawModule from a main source file with preprocessing.
        
        Args:
            main_file: Path to the main source file (which may include other files)
            options: Compiler options for CuPy RawModule
            
        Returns:
            CuPy RawModule and the combined source code
        """
        # Reset line mapping
        self.line_map = {}
        
        # Process the main file (which will recursively process includes)
        processed_source = self.preprocess_file(main_file)
        
        # Create the module
        try:
            if not DEBUG:
                options += ("-DNDEBUG",)
                
            mod = cp.RawModule(code=processed_source, options = options)
            
            # Return a wrapped module that handles function lookup errors
            return self.WrappedCudaModule(mod, processed_source, self), processed_source
        except Exception as e:
            # Process the error message to map line numbers to original files
            error_msg = str(e)
            self.process_error_message(error_msg, processed_source)
            raise
    

    class WrappedCudaModule:
        """Wrapper for CuPy RawModule that provides better error handling."""
        
        def __init__(self, module, source, manager):
            self.module = module
            self.source = source
            self.manager = manager
            
        def get_function(self, name, sizes,kernel_type):
            """Get a kernel function with error handling."""
            try:
                func = self.module.get_function(name)

                return ckw.CudaKernelWrapper(func, sizes, kernel_type=kernel_type)

            except Exception as e:
                print(f"\n{Colors.BOLD}===== CUDA FUNCTION ERROR ====={Colors.RESET}")
                print(f"Error getting function '{Colors.CYAN}{name}{Colors.RESET}':")
                
                error_msg = str(e)
                
                # Process compilation errors in temporary files
                # Pattern matches: /tmp/path/file.cu(line): error message
                # Choose pattern based on platform
                if sys.platform.startswith('win'):  # Windows
                    temp_file_pattern = r'([A-Za-z]:\\.*?\\[^(]+)\((\d+)\):\s*(.+)'
                elif sys.platform.startswith('linux'):  # Linux
                    temp_file_pattern = r'(/tmp/[^(]+)\((\d+)\):\s*(.+)'
                matches = re.finditer(temp_file_pattern, error_msg)
                
                has_line_errors = False
                for match in matches:
                    has_line_errors = True
                    temp_file = match.group(1)
                    line_num = int(match.group(2))
                    error_text = match.group(3)
                    
                    # Determine if this is a warning or error
                    is_warning = "warning" in error_text.lower()
                    is_error = "error" in error_text.lower()
                    message_type = f"{Colors.YELLOW}Warning{Colors.RESET}" if is_warning else f"{Colors.RED}Error{Colors.RESET}"
                    color_prefix = Colors.YELLOW if is_warning else Colors.RED
                    
                    # Map the line number from the temporary file to our combined source
                    if line_num in self.manager.line_map:
                        orig_file, orig_line, include_file = self.manager.line_map[line_num]
                        
                        # For .cu files, adjust the line number by 1 to match the actual error location
                        if orig_file.endswith('.cu'):
                            orig_line = max(1, orig_line - 1)
                        
                        # Get the error line from combined source
                        combined_lines = self.source.split('\n')
                        if 0 <= line_num-1 < len(combined_lines):
                            error_line = combined_lines[line_num-1]
                        else:
                            error_line = "<line not available>"
                        
                        # Print detailed error information
                        if include_file:
                            print(f"{message_type} in file included from {Colors.BLUE}{orig_file}:{Colors.BOLD}{orig_line}{Colors.RESET}")
                            print(f"  Include file: {Colors.BLUE}{include_file}{Colors.RESET}")
                        else:
                            print(f"{message_type} in {Colors.BLUE}{orig_file}:{Colors.BOLD}{orig_line}{Colors.RESET}")
                        
                        print(f"  {Colors.CYAN}{error_line.strip()}{Colors.RESET}")
                        print(f"  {message_type}: {color_prefix}{error_text.strip()}{Colors.RESET}")
                        print()
                
                pass
                # If no line errors were processed, print the raw error
                if not has_line_errors:
                    print(f"  {Colors.RED}{error_msg}{Colors.RESET}")
                
                # Check if it's a missing function error
                if "not found in module" in error_msg:
                    # Search for similar function names to help with typos
                    ptx = self.module.get_info()
                    function_pattern = r'\.visible \.entry ([^\(]+)'
                    available_functions = re.findall(function_pattern, ptx)
                    
                    if available_functions:
                        print(f"\n{Colors.GREEN}Available functions in the module:{Colors.RESET}")
                        for func in sorted(available_functions):
                            print(f"  {Colors.CYAN}{func}{Colors.RESET}")
                    else:
                        print(f"\n{Colors.YELLOW}No functions found in the module. Check for compilation errors.{Colors.RESET}")
                
                print(f"{Colors.BOLD}===== END ERROR ====={Colors.RESET}\n")
                raise
        
        def __getattr__(self, name):
            """Pass through other attributes to the underlying module."""
            return getattr(self.module, name)
    
    def process_error_message(self, error_msg, combined_source):
        """Process CUDA compilation error messages to show original file locations."""
        # Regular expression to find line numbers in error messages
        line_pattern = r'<source>:(\d+):'
        
        # Find all line number references
        matches = re.finditer(line_pattern, error_msg)
        
        print("\n===== CUDA COMPILATION ERROR =====")
        
        # Process each match
        for match in matches:
            line_num = int(match.group(1))
            
            # Get original file and line
            if line_num in self.line_map:
                orig_file, orig_line, include_file = self.line_map[line_num]
                
                # Get the error line from combined source
                combined_lines = combined_source.split('\n')
                if 0 <= line_num-1 < len(combined_lines):
                    error_line = combined_lines[line_num-1]
                else:
                    error_line = "<line not available>"
                
                # Print detailed error information
                if include_file:
                    print(f"Error in file included from {orig_file}:{orig_line}")
                    print(f"  Include file: {include_file}")
                else:
                    print(f"Error in {orig_file}:{orig_line}")
                
                print(f"  {error_line.strip()}")
                
                # Extract the specific error message
                error_line_match = re.search(f'<source>:{line_num}:(.*?)(?:\n|$)', error_msg)
                if error_line_match:
                    print(f"  Error: {error_line_match.group(1).strip()}")
                
                print()
        
        # Print the original error message
        print("Original error message:")
        print(error_msg)
        print("===== END ERROR =====\n")

    def save_combined_source(self, combined_source, output_file):
        """Save the combined source to a file for debugging."""
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(combined_source)
            
        # Also save a version with line numbers for easier debugging
        with open(f"{output_file}.numbered", 'w', encoding='utf-8') as f:
            for i, line in enumerate(combined_source.split('\n'), 1):
                f.write(f"{i:5d}: {line}\n") 