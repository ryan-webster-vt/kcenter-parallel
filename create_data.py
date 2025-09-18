import random

# Parameters
num_clusters = 4       # number of clusters
points_per_cluster = 50  # number of points in each cluster
spread = 1.0           # standard deviation around cluster center
output_file = "data.txt"

# Generate cluster centers randomly
cluster_centers = [(random.uniform(0, 20), random.uniform(0, 20)) for _ in range(num_clusters)]

points = []

for cx, cy in cluster_centers:
    for _ in range(points_per_cluster):
        x = random.gauss(cx, spread)
        y = random.gauss(cy, spread)
        points.append((x, y))

# Save to file
with open(output_file, "w") as f:
    for x, y in points:
        f.write(f"{x:.4f} {y:.4f}\n")

print(f"Generated {num_clusters * points_per_cluster} points in '{output_file}'")
