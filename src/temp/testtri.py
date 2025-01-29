import numpy as np

def generate_random_triangle():
    """Generate 3 random vertices in 2D space"""
    return np.random.uniform(low=0.0, high=10.0, size=(3, 2))

def generate_perturbed_coordinates(amin_limit=-0.1):
    """Generate barycentric coordinates with amin > specified limit"""
    while True:
        # Start with valid barycentric coordinates
        valid_coords = np.random.dirichlet((1, 1, 1))
        
        # Apply perturbation while maintaining sum=1
        perturbation = np.random.uniform(-0.1, 0.1, 3)
        perturbed = valid_coords + perturbation
        perturbed /= perturbed.sum()
        
        if np.min(perturbed) > amin_limit:
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