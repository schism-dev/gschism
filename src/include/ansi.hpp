// ansi.hpp
namespace ansi {
    // Colors
    std::string RESET = "\033[0m";
    std::string RED = "\033[31m";
    std::string GREEN = "\033[32m";
    std::string LIGHT_GREEN = "\033[92m";
    std::string BLUE = "\033[34m";
    std::string LIGHT_BLUE = "\033[94m";
    std::string YELLOW = "\033[33m";
    std::string WHITE = "\033[37m";
    std::string CYAN = "\033[36m";
    std::string MAGENTA = "\033[35m";
    std::string BLACK = "\033[30m";
    std::string GRAY = "\033[90m";
    std::string LIGHT_GRAY = "\033[37m";
    std::string BOLD = "\033[1m";

    inline std::string red(const std::string& text) { return RED + text + RESET; }
    inline std::string green(const std::string& text) { return GREEN + text + RESET; }
    inline std::string blue(const std::string& text) { return BLUE + text + RESET; }
    inline std::string yellow(const std::string& text) { return YELLOW + text + RESET; }
    inline std::string white(const std::string& text) { return WHITE + text + RESET; }
    inline std::string cyan(const std::string& text) { return CYAN + text + RESET; }
    inline std::string magenta(const std::string& text) { return MAGENTA + text + RESET; }
    inline std::string black(const std::string& text) { return BLACK + text + RESET; }
    inline std::string gray(const std::string& text) { return GRAY + text + RESET; }
    inline std::string light_gray(const std::string& text) { return LIGHT_GRAY + text + RESET; }  
    inline std::string light_green(const std::string& text) { return LIGHT_GREEN + text + RESET; }
    inline std::string light_blue(const std::string& text) { return LIGHT_BLUE + text + RESET; }

    inline std::string bold(const std::string& text) { return BOLD + text + RESET; }
    // Cursor movement
    inline std::string cursorUp(int n = 1) { return "\033[" + std::to_string(n) + "A"; }
    inline std::string cursorDown(int n = 1) { return "\033[" + std::to_string(n) + "B"; }
    inline std::string cursorForward(int n = 1) { return "\033[" + std::to_string(n) + "C"; }
    inline std::string cursorBack(int n = 1) { return "\033[" + std::to_string(n) + "D"; }
    inline std::string cursorTo(int row, int col) { return "\033[" + std::to_string(row) + ";" + std::to_string(col) + "H"; }
    
    // Clear screen operations
    inline std::string clearScreen() { return "\033[2J"; }
    inline std::string clearLine() { return "\033[2K"; }
    inline std::string saveCursor() { return "\033[s"; }
    inline std::string restoreCursor() { return "\033[u"; }

        // Unicode block characters for sub-character precision
    const char* const blocks[] = {" ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"};

    inline std::string progressBar(int current, int total, int barWidth = 20) {
        float progress = static_cast<float>(current) / total;
        float fullBlocks = barWidth * progress;
        int numFullBlocks = static_cast<int>(fullBlocks);
        int subBlock = static_cast<int>((fullBlocks - numFullBlocks) * 8);
        
        std::string bar = "[";
        
        // Add full blocks
        for (int i = 0; i < numFullBlocks; ++i) {
            bar += blocks[8];  // Full block
        }
        
        // Add sub-block if any
        if (numFullBlocks < barWidth) {
            bar += blocks[subBlock];
            // Fill the rest with spaces
            for (int i = numFullBlocks + 1; i < barWidth; ++i) {
                bar += blocks[0];  // Space
            }
        }
        
        int percent = static_cast<int>(progress * 100.0f);
        bar += "] " + std::to_string(percent) + "% ";
        bar += "(" + std::to_string(current) + "/" + std::to_string(total) + ")";
        
        return bar;
    }

    // Colored version
    inline std::string progressBar(int current, int total, const std::string& color, int barWidth = 20) {
        return color + progressBar(current, total, barWidth) + RESET;
    }

    // Version with description
    inline std::string progressBar(const std::string& desc, int current, int total, 
                                 const std::string& color, int barWidth = 20) {
        std::string paddedDesc = desc;
        if (paddedDesc.length() < 45) {
            paddedDesc.insert(0, 45 - paddedDesc.length(), ' ');
        }
        return color + paddedDesc + ": " + progressBar(current, total, barWidth) + RESET;
    }    
}


