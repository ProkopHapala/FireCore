
import os
import sys
import csv
import json
import glob

def find_latest_scan_dir(base_dir="scan_outputs"):
    if not os.path.exists(base_dir):
        return None
    subdirs = [os.path.join(base_dir, d) for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    if not subdirs:
        return None
    return max(subdirs, key=os.path.getmtime)

def generate_html(scan_dir):
    csv_path = os.path.join(scan_dir, "summary.csv")
    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} not found.")
        return

    data = []
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Convert numeric fields
            for key in ['bmix', 'dt', 'relax', 'alpha', 'nsteps', 'it_1e2', 'it_1e3', 'last_max_err']:
                try:
                    if key in row:
                        row[key] = float(row[key]) if '.' in row[key] or 'e' in row[key] else int(row[key])
                except ValueError:
                    pass # Keep as string if conversion fails (e.g. status)
            data.append(row)
            
    # Format data as a JS string
    json_data = json.dumps(data, indent=2)

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>XPDB Convergence Analysis</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/d3/7.9.0/d3.js"></script>
    <style>
        body {{ font-family: sans-serif; margin: 20px; }}
        .controls {{ margin-bottom: 20px; padding: 10px; background: #f0f0f0; border-radius: 5px; }}
        .chart-container {{ display: flex; flex-wrap: wrap; }}
        .chart {{ margin: 10px; border: 1px solid #ddd; padding: 10px; }}
        .tooltip {{ position: absolute; background: white; border: 1px solid #333; padding: 5px; pointer-events: none; opacity: 0; transition: opacity 0.2s; }}
        label {{ margin-right: 15px; }}
    </style>
</head>
<body>
    <h1>XPDB Convergence Analysis</h1>
    <div class="controls">
        <label><strong>Y-Axis Metric:</strong> 
            <select id="metricSelect">
                <option value="it_1e2" selected>Iterations to 1e-2</option>
                <option value="it_1e3">Iterations to 1e-3</option>
                <option value="it_dx_1e2">Iters Translation 1e-2</option>
                <option value="it_dx_1e3">Iters Translation 1e-3</option>
                <option value="it_dang_1e2">Iters Rotation 1e-2</option>
                <option value="it_dang_1e3">Iters Rotation 1e-3</option>
                <option value="last_max_err">Final Max Error</option>
                <option value="last_max_dx">Final Max dx</option>
                <option value="last_max_dang">Final Max dang</option>
            </select>
        </label>
        <label><strong>Relaxation Filter:</strong>
            <select id="relaxSelect">
                <option value="all">All</option>
            </select>
        </label>
        <label><strong>Kernel Filter:</strong>
            <select id="kernelSelect">
                <option value="all">All</option>
            </select>
        </label>
    </div>
    
    <div id="charts"></div>
    <div class="tooltip" id="tooltip"></div>

    <script>
        // Embedded Data
        const rawData = {json_data};

        // Parse and prepare data
        rawData.forEach(d => {{
            d.bmix = +d.bmix;
            d.dt = +d.dt;
            d.relax = +d.relax;
            d.it_1e2 = +d.it_1e2;
            d.it_1e3 = +d.it_1e3;
            d.it_dx_1e2 = +d.it_dx_1e2 || -1;
            d.it_dx_1e3 = +d.it_dx_1e3 || -1;
            d.it_dang_1e2 = +d.it_dang_1e2 || -1;
            d.it_dang_1e3 = +d.it_dang_1e3 || -1;
            
            // Handle non-convergence (-1)
            d.it_1e2_plot = d.it_1e2 < 0 ? 600 : d.it_1e2; // Cap at 600 for plot visibility
            d.it_1e3_plot = d.it_1e3 < 0 ? 600 : d.it_1e3;
            d.it_dx_1e2_plot = d.it_dx_1e2 < 0 ? 600 : d.it_dx_1e2;
            d.it_dx_1e3_plot = d.it_dx_1e3 < 0 ? 600 : d.it_dx_1e3;
            d.it_dang_1e2_plot = d.it_dang_1e2 < 0 ? 600 : d.it_dang_1e2;
            d.it_dang_1e3_plot = d.it_dang_1e3 < 0 ? 600 : d.it_dang_1e3;
        }});

        // Get unique values
        const relaxValues = [...new Set(rawData.map(d => d.relax))].sort((a,b) => a-b);
        const dtValues = [...new Set(rawData.map(d => d.dt))].sort((a,b) => a-b);
        const kernelValues = [...new Set(rawData.map(d => d.kernel))].sort();

        // Populate controls
        const relaxSelect = d3.select("#relaxSelect");
        relaxValues.forEach(r => {{
            relaxSelect.append("option").text(r).attr("value", r);
        }});
        // Select 1.0 by default if available
        if(relaxValues.includes(1.0)) relaxSelect.property("value", "1.0");

        const kernelSelect = d3.select("#kernelSelect");
        kernelValues.forEach(k => {{
            kernelSelect.append("option").text(k).attr("value", k);
        }});

        function drawCharts() {{
            const metric = d3.select("#metricSelect").property("value");
            const relaxFilter = d3.select("#relaxSelect").property("value");
            const kernelFilter = d3.select("#kernelSelect").property("value");
            
            d3.select("#charts").html(""); // Clear existing

            // Filter data
            let plotData = rawData.filter(d => d.status === 'ok');
            if (relaxFilter !== "all") {{
                plotData = plotData.filter(d => d.relax == relaxFilter);
            }}
            if (kernelFilter !== "all") {{
                plotData = plotData.filter(d => d.kernel == kernelFilter);
            }}

            // If "all" relaxation, maybe group by dt?
            // Let's make one main chart: X=bmix, Y=metric, Series=dt
            
            const margin = {{top: 20, right: 80, bottom: 50, left: 60}};
            const width = 800 - margin.left - margin.right;
            const height = 500 - margin.top - margin.bottom;

            const svg = d3.select("#charts").append("svg")
                .attr("width", width + margin.left + margin.right)
                .attr("height", height + margin.top + margin.bottom)
                .append("g")
                .attr("transform", `translate(${{margin.left}},${{margin.top}})`);

            // Scales
            const xExtent = d3.extent(plotData, d => d.bmix);
            const xScale = d3.scaleLinear()
                .domain([Math.max(0, xExtent[0] - 0.05), Math.min(1, xExtent[1] + 0.05)])
                .range([0, width]);

            let yDomain;
            if (metric.startsWith("it")) {{
                 yDomain = [0, 550]; // Fixed scale for iterations
            }} else {{
                 yDomain = [0, d3.max(plotData, d => d[metric]) * 1.1];
            }}
            
            const yScale = d3.scaleLinear()
                .domain(yDomain)
                .range([height, 0]);

            const colorScale = d3.scaleOrdinal()
                .domain(dtValues)
                .range(d3.schemeCategory10);

            // Axes
            svg.append("g")
                .attr("transform", `translate(0,${{height}})`)
                .call(d3.axisBottom(xScale).ticks(10));
            
            svg.append("text")
                .attr("x", width/2)
                .attr("y", height + 40)
                .style("text-anchor", "middle")
                .text("bmix (Momentum)");

            svg.append("g")
                .call(d3.axisLeft(yScale));

            svg.append("text")
                .attr("transform", "rotate(-90)")
                .attr("y", -45)
                .attr("x", -height/2)
                .style("text-anchor", "middle")
                .text(metric === 'last_max_err' ? 'Max Error [A]' : 'Iterations');

            // Group by dt to draw lines
            const grouped = d3.group(plotData, d => d.dt);
            
            // Draw lines
            grouped.forEach((values, dt) => {{
                // Sort by bmix for line connection
                values.sort((a,b) => a.bmix - b.bmix);
                
                const line = d3.line()
                    .x(d => xScale(d.bmix))
                    .y(d => metric.startsWith("it") ? yScale(d[metric + "_plot"]) : yScale(d[metric]));

                svg.append("path")
                    .datum(values)
                    .attr("fill", "none")
                    .attr("stroke", colorScale(dt))
                    .attr("stroke-width", 2)
                    .attr("d", line);
            }});

            // Draw points
            svg.selectAll(".dot")
                .data(plotData)
                .enter().append("circle")
                .attr("class", "dot")
                .attr("cx", d => xScale(d.bmix))
                .attr("cy", d => metric.startsWith("it") ? yScale(d[metric + "_plot"]) : yScale(d[metric]))
                .attr("r", 5)
                .attr("fill", d => colorScale(d.dt))
                .on("mouseover", function(event, d) {{
                    d3.select("#tooltip")
                        .style("opacity", 1)
                        .html(`
                            <strong>${{d.tag}}</strong><br>
                            bmix: ${{d.bmix}}<br>
                            dt: ${{d.dt}}<br>
                            relax: ${{d.relax}}<br>
                            alpha: ${{d.alpha.toExponential(2)}}<br>
                            ${{metric}}: ${{d[metric]}}
                        `)
                        .style("left", (event.pageX + 10) + "px")
                        .style("top", (event.pageY - 28) + "px");
                }})
                .on("mouseout", function() {{
                    d3.select("#tooltip").style("opacity", 0);
                }});

            // Legend
            const legend = svg.append("g")
                .attr("transform", `translate(${{width - 60}}, 0)`);
            
            dtValues.forEach((dt, i) => {{
                legend.append("rect")
                    .attr("x", 0)
                    .attr("y", i * 20)
                    .attr("width", 10)
                    .attr("height", 10)
                    .style("fill", colorScale(dt));
                
                legend.append("text")
                    .attr("x", 15)
                    .attr("y", i * 20 + 9)
                    .text(`dt=${{dt}}`)
                    .style("font-size", "12px")
                    .attr("alignment-baseline", "middle");
            }});
        }}

        // Event listeners
        d3.select("#metricSelect").on("change", drawCharts);
        d3.select("#relaxSelect").on("change", drawCharts);
        d3.select("#kernelSelect").on("change", drawCharts);

        // Initial draw
        drawCharts();

    </script>
</body>
</html>"""

    out_path = os.path.join(scan_dir, "convergence_viz.html")
    with open(out_path, 'w') as f:
        f.write(html_content)
    
    print(f"Visualization generated: {out_path}")

if __name__ == "__main__":
    scan_dir = None
    if len(sys.argv) > 1:
        scan_dir = sys.argv[1]
    else:
        scan_dir = find_latest_scan_dir()
    
    if scan_dir:
        print(f"Processing scan directory: {scan_dir}")
        generate_html(scan_dir)
    else:
        print("No scan directory found.")
