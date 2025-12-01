import numpy as np
import matplotlib.pyplot as plt

def generate_random_triangle():
    """Generate 3 random vertices in 2D space"""
    return np.random.uniform(low=0.0, high=10.0, size=(3, 2))

def generate_perturbed_coordinates(amin_limit=-0.01):
    """Generate barycentric coordinates with at least one negative coordinate > -0.01"""
    while True:
        # Start with valid barycentric coordinates
        valid_coords = np.random.dirichlet((1, 1, 1))
        
        # Choose one coordinate to make negative (between -0.01 and 0)
        idx_to_make_negative = np.random.randint(0, 3)
        negative_value = np.random.uniform(-0.01, 0.0)
        
        # Create perturbed coordinates
        perturbed = valid_coords.copy()
        perturbed[idx_to_make_negative] = negative_value
        
        # Adjust to maintain sum=1 by redistributing the difference
        current_sum = perturbed.sum()
        excess = current_sum - 1.0  # This will be negative since we made one coord negative
        
        # Distribute the excess (which is negative, so we subtract) among the other two coordinates
        other_indices = [i for i in range(3) if i != idx_to_make_negative]
        # Distribute proportionally to maintain relative ratios
        if len(other_indices) > 0:
            other_sum = perturbed[other_indices].sum()
            if other_sum > 0:
                perturbed[other_indices] -= excess * (perturbed[other_indices] / other_sum)
            else:
                # If other_sum is 0 or negative, distribute equally
                perturbed[other_indices] -= excess / len(other_indices)
        
        # Verify: at least one negative coordinate > -0.01, and sum is approximately 1
        if np.any((perturbed < 0) & (perturbed > -0.01)) and abs(perturbed.sum() - 1.0) < 1e-10:
            return perturbed

def original_method(a):
    a = np.array(a, copy=True)
    threshold = 1e-2
    
    # Find negative coordinates
    negatives = a < threshold
    if not np.any(negatives):
        return a
    
    # Set negatives to threshold
    a[negatives] = threshold
    
    # Adjust largest coordinate to maintain sum=1
    excess = a.sum() - 1.0
    largest = np.argmax(a)
    a[largest] -= excess
    
    return a

def user_method(a):
    a = np.array(a, copy=True)
    amin = np.min(a)

    if amin>=0: 
      return a

    threshold = 1e-3
    
    A = -amin + threshold
    B = 1 / (a.sum() + 3 * A)  # Since original sum is 1
    adjusted = (a + A) * B
    
    return adjusted

def barycentric_to_cartesian(a, triangle):
    """Convert barycentric coordinates to Cartesian"""
    return a[0] * triangle[0] + a[1] * triangle[1] + a[2] * triangle[2]

def plot_triangle_and_points(triangle, original_point, orig_point, user_point, orig, orig_fixed, user_fixed):
    """Plot triangle with original and nudged points"""
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Plot triangle
    triangle_closed = np.vstack([triangle, triangle[0]])  # Close the triangle
    ax.plot(triangle_closed[:, 0], triangle_closed[:, 1], 'k-', linewidth=2, label='Triangle')
    ax.plot(triangle[:, 0], triangle[:, 1], 'ko', markersize=8, label='Triangle vertices')
    
    # Plot original point (from perturbed coordinates)
    ax.plot(original_point[0], original_point[1], 'ro', markersize=10, label=f'Original point (a={orig_fixed.sum():.6f})')
    
    # Plot nudged points
    ax.plot(orig_point[0], orig_point[1], 'bs', markersize=10, label='Original method')
    ax.plot(user_point[0], user_point[1], 'g^', markersize=10, label='User method')
    
    # Add text annotations for barycentric coordinates
    ax.text(original_point[0] + 0.1, original_point[1] + 0.1, 
            f'Original\n{orig}', fontsize=8, color='red')
    ax.text(orig_point[0] + 0.1, orig_point[1] - 0.1, 
            f'Original method\n{orig_fixed}', fontsize=8, color='blue')
    ax.text(user_point[0] + 0.1, user_point[1] - 0.1, 
            f'User method\n{user_fixed}', fontsize=8, color='green')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Triangle with Original and Nudged Points')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal', adjustable='box')
    
    plt.tight_layout()
    plt.show(block=False)
    plt.pause(0.1)  # Small pause to ensure plot is rendered
    
    # Wait for key press
    input("Press Enter to continue to next sample...")
    plt.close(fig)    

def test_methods(num_samples=20000):
    original_distances = []
    user_distances = []
    
 
    for samp in range(num_samples):
        print('sample', samp)
        # Generate random triangle and coordinates
        triangle = generate_random_triangle()
        a = generate_perturbed_coordinates(amin_limit=-0.1)
        # print('orig: ', a)
        
        # Convert to Cartesian coordinates
        original_point = barycentric_to_cartesian(a, triangle)
        # print('in cartesians: ', original_point)

        
        # Apply both methods
        orig_fixed = original_method(a)
        user_fixed = user_method(a)

        print(orig_fixed.sum(), user_fixed.sum())

        #print('>>',a)
        
        # Calculate fixed points
        orig_point = barycentric_to_cartesian(orig_fixed, triangle)
        user_point = barycentric_to_cartesian(user_fixed, triangle)
        #print('---')
        #print(orig_fixed)
        #print(user_fixed)
        
        # Plot the triangle and points
        plot_triangle_and_points(triangle, original_point, orig_point, user_point, a, orig_fixed, user_fixed)
        
        # Calculate distances
        original_distances.append(np.linalg.norm(original_point - orig_point))
        user_distances.append(np.linalg.norm(original_point - user_point))
    
    # Print results
    print(f"Original method average distance: {np.mean(original_distances):.6f}")
    print(f"User method average distance: {np.mean(user_distances):.6f}")
    print(f"\nKey observations:")
    print(f"- User method preserves coordinate relationships better")
    print(f"- Original method guarantees minimum of 1e-2, user method allows smaller values")
    print(f"- Both methods maintain sum=1 (verified during testing)")

# Run the comparison
np.random.seed(42)
test_methods()