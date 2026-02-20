import re

with open("test_assembly.py", "r") as f:
    content = f.read()

# Replace the export mask filtering block
old_block_1 = """    export_mask = (z_spans < args.zspan_max) & (scores < args.clash_max)
    export_indices = np.where(export_mask)[0]
    
    print(f"Configurations passing export criteria (zspan < {args.zspan_max}, clash < {args.clash_max}): {len(export_indices)}")
    
    if export_indices.size > 0:
        export_sorted = export_indices[np.argsort(total_scores[export_indices])]
        export_sorted = export_sorted[:args.export_max] # Cap export size
        best_indices = export_sorted[:args.top_k]
    elif len(valid_scores) == 0:
        best_indices = [int(np.argmin(total_scores))]
        export_sorted = np.array(best_indices, dtype=int)
    else:
        valid_indices = np.where(valid_mask)[0]
        sorted_valid = valid_indices[np.argsort(total_scores[valid_indices])]
        best_indices = sorted_valid[:args.top_k]
        export_sorted = sorted_valid[:args.export_max]
        
    best_idx = best_indices[0]
    print(f"Best configuration index: {best_idx} with Clash {scores[best_idx]:.4f}, Z-span {z_spans[best_idx]:.4f}, Total Score {total_scores[best_idx]:.4f}")"""

new_block_1 = """    export_mask = (z_spans < args.zspan_max) & (scores < args.clash_max)
    export_indices = np.where(export_mask)[0]
    
    print(f"Configurations passing export criteria (zspan < {args.zspan_max}, clash < {args.clash_max}): {len(export_indices)}")
    
    if export_indices.size == 0:
        print(f"WARNING: No configurations passed the strict export criteria (zspan < {args.zspan_max}, clash < {args.clash_max}).")
        print("         The Pareto plot will still be generated, but no XYZ movie or 2D plots will be exported.")
        export_sorted = np.array([], dtype=int)
        best_indices = np.array([], dtype=int)
        best_idx = int(np.argmin(total_scores)) # best overall, even if invalid, for plot reference
    else:
        export_sorted = export_indices[np.argsort(total_scores[export_indices])]
        export_sorted = export_sorted[:args.export_max] # Cap export size
        best_indices = export_sorted[:args.top_k]
        best_idx = best_indices[0]
        
    print(f"Best configuration overall: index {best_idx} with Clash {scores[best_idx]:.4f}, Z-span {z_spans[best_idx]:.4f}, Total Score {total_scores[best_idx]:.4f}")"""

# Replace the Pareto plot code block
old_block_2 = """        # Pareto Plot
        plt.figure(figsize=(8, 6))
        plot_mask = export_mask
        sc = plt.scatter(scores[plot_mask], z_spans[plot_mask], c=total_scores[plot_mask], cmap='viridis', alpha=0.6, s=10)
        plt.colorbar(sc, label=f'Total Score (clash + {args.zpenalty}*zspan)')
        plt.xlabel('Clash Penalty')
        plt.ylabel('Z-Span (A)')
        plt.title('Pareto Frontier: Clash vs Z-Span')
        plt.scatter(scores[best_idx], z_spans[best_idx], color='red', marker='x', s=150, linewidths=3, label='Best Config')
        plt.xlim(-1, min(args.penalty, float(np.max(scores[plot_mask]) * 1.1) if np.any(plot_mask) else args.penalty))
        plt.legend()
        pareto_name = os.path.join(args.outdir, f'pareto_{args.preset}.png')
        plt.savefig(pareto_name)
        print(f"Saved Pareto plot to {pareto_name}")
        plt.close()"""

new_block_2 = """        # Pareto Plot
        plt.figure(figsize=(8, 6))
        
        # Draw all points as small black dots
        plt.scatter(scores, z_spans, c='black', alpha=0.3, s=5, label='All Configs')
        
        # Highlight exported configurations as slightly larger red dots
        if len(export_sorted) > 0:
            plt.scatter(scores[export_sorted], z_spans[export_sorted], c='red', alpha=0.8, s=20, label='Exported')
            plt.scatter(scores[best_idx], z_spans[best_idx], color='lime', marker='*', s=150, linewidths=2, label='Best Exported')
        else:
            plt.scatter(scores[best_idx], z_spans[best_idx], color='orange', marker='x', s=150, linewidths=2, label='Best Overall (Failing)')
            
        # Compute and plot Pareto front
        pts = np.vstack((scores, z_spans)).T
        is_pareto = np.ones(pts.shape[0], dtype=bool)
        for i, p in enumerate(pts):
            if is_pareto[i]:
                dominating = np.all(pts <= p, axis=1) & np.any(pts < p, axis=1)
                is_pareto[i] = not np.any(dominating)
                
        pareto_idx = np.where(is_pareto)[0]
        # Sort pareto front by clash score for line plotting
        pareto_sorted = pareto_idx[np.argsort(scores[pareto_idx])]
        plt.plot(scores[pareto_sorted], z_spans[pareto_sorted], 'b-', alpha=0.5, lw=2, label='Pareto Front')
        
        # Plot thresholds if within limits
        plt.axvline(x=args.clash_max, color='r', linestyle='--', alpha=0.3, label='Clash Max')
        plt.axhline(y=args.zspan_max, color='b', linestyle='--', alpha=0.3, label='Z-Span Max')
        
        plt.xlabel('Clash Penalty')
        plt.ylabel('Z-Span (A)')
        plt.title('Pareto Frontier: Clash vs Z-Span')
        
        # Set limits dynamically but reasonably
        max_clash_plot = max(args.clash_max * 2, np.percentile(scores, 10)) # Focus on the good 10% or at least twice the max
        max_zspan_plot = max(args.zspan_max * 1.2, np.percentile(z_spans, 50))
        plt.xlim(-0.1, min(args.penalty, max_clash_plot))
        plt.ylim(0, max_zspan_plot)
        
        plt.legend()
        pareto_name = os.path.join(args.outdir, f'pareto_{args.preset}.png')
        plt.savefig(pareto_name, dpi=150)
        print(f"Saved Pareto plot to {pareto_name}")
        plt.close()"""

content = content.replace(old_block_1, new_block_1)
content = content.replace(old_block_2, new_block_2)

# Fix indentation of dump blocks
old_dump_loop = """        # Dump configuration movie
        out_name = os.path.join(args.outdir, f"assembly_best_{args.preset}_movie.xyz")
        print(f"Dumping {len(export_sorted)} configurations (capped at {args.export_max}) to {out_name} with {nmols_out} replicas...")
        with open(out_name, 'w') as f:
            for i_rank, idx in enumerate(export_sorted):
                best_transforms = transforms[idx][:nmols_out]"""

new_dump_loop = """        # Dump configuration movie
        if len(export_sorted) > 0:
            out_name = os.path.join(args.outdir, f"assembly_best_{args.preset}_movie.xyz")
            print(f"Dumping {len(export_sorted)} configurations (capped at {args.export_max}) to {out_name} with {nmols_out} replicas...")
            with open(out_name, 'w') as f:
                for i_rank, idx in enumerate(export_sorted):
                    best_transforms = transforms[idx][:nmols_out]"""

content = content.replace(old_dump_loop, new_dump_loop)

# Fix indentation of out_atoms loop manually
old_out_atoms = """                for i in range(len(out_atoms)):
                    orig_idx = i % mol.natoms
                    elem = mol.enames[orig_idx]
                    if args.xyz_simple:
                        elem = elem.split('_')[0] if '_' in elem else elem
                    f.write(f"{elem} {out_atoms[i,0]:.6f} {out_atoms[i,1]:.6f} {out_atoms[i,2]:.6f}\\n")
        print(f"Saved movie to {out_name}")"""

new_out_atoms = """                for i in range(len(out_atoms)):
                    orig_idx = i % mol.natoms
                    elem = mol.enames[orig_idx]
                    if args.xyz_simple:
                        elem = elem.split('_')[0] if '_' in elem else elem
                    f.write(f"{elem} {out_atoms[i,0]:.6f} {out_atoms[i,1]:.6f} {out_atoms[i,2]:.6f}\\n")
            print(f"Saved movie to {out_name}")"""
            
content = content.replace(old_out_atoms, new_out_atoms)

old_plot_loop = """        # Render Atoms Plot for the top K configs
        plot_k = min(args.plot_best_k, len(best_indices))
        for rank in range(plot_k):
            idx = best_indices[rank]"""

new_plot_loop = """        # Render Atoms Plot for the top K configs
        plot_k = min(args.plot_best_k, len(best_indices))
        for rank in range(plot_k):
            idx = best_indices[rank]"""

content = content.replace(old_plot_loop, new_plot_loop)


with open("test_assembly.py", "w") as f:
    f.write(content)

