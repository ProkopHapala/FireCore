import argparse
from ase.io import read, write
from mace.calculators import MACECalculator
import matplotlib.pyplot as plt
import numpy as np
import faiss
from tqdm import tqdm
import os

def compute_descriptors(structures, calculator):
    descriptors = []
    for i, atoms in enumerate(tqdm(structures, desc="Computing descriptors")):
        try:
            desc = calculator.get_descriptors(atoms, invariants_only=False)
            desc_mean = np.mean(desc, axis=0)  # Mean-pool to fixed-size descriptor
            descriptors.append(desc_mean)
        except Exception as e:
            print(f"Error at structure {i}: {e}")
    return np.vstack(descriptors).astype(np.float32)

def filter_descriptors(descriptor_matrix, threshold):
    D = descriptor_matrix.shape[1]
    index = faiss.IndexFlatL2(D)
    accepted_indices = []

    for i in tqdm(range(len(descriptor_matrix)), desc="Filtering with FAISS"):
        desc = descriptor_matrix[i].reshape(1, -1)
        if index.ntotal == 0:
            index.add(desc)
            accepted_indices.append(i)
            continue
        dist, _ = index.search(desc, 1)
        if dist[0][0] >= threshold**2:
            index.add(desc)
            accepted_indices.append(i)
    return accepted_indices

def sanity_check_nn(descriptors, structures, query_index, k, output_file="nearest_neighbors.extxyz"):
    print(f"\nRunning sanity check: finding {k} nearest neighbors for structure {query_index}")
    index = faiss.IndexFlatL2(descriptors.shape[1])
    index.add(descriptors)

    query_desc = descriptors[query_index].reshape(1, -1)
    distances, indices = index.search(query_desc, k)

    print(f"🔍 Nearest neighbors (index, distance):")
    for dist, idx in zip(distances[0], indices[0]):
        print(f"  - Index: {idx}, Distance: {dist:.4f}")

    to_save = [structures[i] for i in indices[0]]
    write(output_file, to_save)
    print(f"Nearest neighbors written to: {output_file}")

def plot_distance_histogram(descriptors, output_path="distance_histogram.png", cdf_output_path="distance_cdf.png"):
    print("📊 Computing nearest neighbor distances for histogram and CDF...")
    index = faiss.IndexFlatL2(descriptors.shape[1])
    index.add(descriptors)

    # For each descriptor, get the distance to its nearest other structure (excluding itself)
    distances, _ = index.search(descriptors, 2)
    nearest_dists = np.sqrt(distances[:, 1])  # Skip self-match

    # Histogram plot
    print(f"📈 Plotting histogram of {len(nearest_dists)} distances...")
    plt.figure(figsize=(8, 5))
    plt.hist(nearest_dists, bins=100, color='skyblue', edgecolor='black')
    plt.xlabel("L2 distance to nearest neighbor")
    plt.ylabel("Frequency")
    plt.title("Descriptor Space: Nearest Neighbor Distances")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Distance histogram saved to: {output_path}")

    # CDF plot
    sorted_dists = np.sort(nearest_dists)
    cumulative = np.arange(1, len(sorted_dists)+1) / len(sorted_dists)

    plt.figure(figsize=(8, 5))
    plt.plot(sorted_dists, cumulative, color='darkorange')
    plt.xlabel("L2 distance to nearest neighbor")
    plt.ylabel("Cumulative fraction")
    plt.title("CDF of Nearest Neighbor Distances")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(cdf_output_path)
    print(f"CDF plot saved to: {cdf_output_path}")

def main():
    parser = argparse.ArgumentParser(description="Filter .extxyz dataset using MACE descriptors and FAISS")
    parser.add_argument("--input", required=True, help="Input .extxyz file with structures")
    parser.add_argument("--model", required=True, help="Path to MACE model (.model)")
    parser.add_argument("--output", default="filtered.extxyz", help="Output file for filtered structures")
    parser.add_argument("--threshold", type=float, default=1.0, help="Descriptor distance threshold")
    parser.add_argument("--device", default="cuda", help="Device to run MACE on (cuda or cpu)")
    parser.add_argument("--batch_size", type=int, default=256, help="Batch size for descriptor calculation")
    parser.add_argument("--plot_distances", action="store_true", help="Plot histogram of nearest neighbor distances")
    parser.add_argument("--hist_output", default="distance_histogram.png", help="Output path for distance histogram image")

    # Nearest neighbor check
    parser.add_argument("--nn_index", type=int, help="Index of structure for nearest neighbor sanity check")
    parser.add_argument("--nn_k", type=int, default=5, help="Number of neighbors to retrieve")
    parser.add_argument("--nn_output", default="nearest_neighbors.extxyz", help="Output file for nearest neighbors")

    args = parser.parse_args()

    print(f"Reading structures from {args.input}")
    structures = read(args.input, ":")

    print(f"Loading MACE model from {args.model}")
    calculator = MACECalculator(model_paths=args.model, device=args.device)

    print("Computing descriptors...")
    descriptors = compute_descriptors(structures, calculator)

    if args.plot_distances:
        cdf_path = os.path.splitext(args.hist_output)[0] + "_cdf.png"
        plot_distance_histogram(descriptors, args.hist_output, cdf_path)

    if args.nn_index is not None:
        sanity_check_nn(descriptors, structures, args.nn_index, args.nn_k, args.nn_output)

    print("Filtering structures with FAISS...")
    accepted_indices = filter_descriptors(descriptors, threshold=args.threshold)

    print(f"{len(accepted_indices)} out of {len(structures)} structures accepted.")
    filtered_structures = [structures[i] for i in accepted_indices]
    write(args.output, filtered_structures)
    print(f"Filtered structures written to {args.output}")

if __name__ == "__main__":
    main()

