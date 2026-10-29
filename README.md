# Colombo Road Network: A Graph-Theoretic Study of Flood Resilience

**Project focus:** Graph Theory · Network Analysis · Urban Road Networks · Infrastructure Resilience

## Overview

This project applies graph theory to study the structure and resilience of the Colombo road network. Using road-network data from OpenStreetMap, the analysis examines network connectivity, bridge edges, edge betweenness centrality, and budget-constrained road-link addition scenarios.

The objective is to understand how network structure can help identify vulnerable connections and explore strategies for improving connectivity.

> **Important limitation:** This is a structural network analysis, not a flood-damage prediction model. Flood-event data and validated flood-impact layers were not available for this analysis.

## Objectives

- Represent the Colombo road network as a graph.
- Identify bridge edges and examine their effect on network connectivity.
- Analyze edge betweenness centrality to identify structurally important connections.
- Explore alternative link-addition scenarios under budget constraints.
- Compare network robustness and connectivity-related outcomes.

## Technologies Used

- **Python** — data analysis and implementation
- **NetworkX** — graph algorithms and network analysis
- **Pandas & NumPy** — data processing
- **Matplotlib** — visualization
- **SciPy** — scientific computing and optimization support
- **OpenStreetMap** — road-network data source
- **Google Colab / Jupyter Notebook** — exploratory analysis

## Repository Structure

```text
colombo-flood-resilience-graph-theory/
├── README.md
├── colombo_analysis.py
├── colombo_nodes.csv
├── colombo_roads.csv
├── notebook_original_colab.ipynb
├── fig1_bridge_map.png
├── fig2_impact_betweenness.png
├── fig3_edi.png
├── fig4_scenarios.png
├── fig5_augmentation_map.png
├── fig6_budget_curve.png
└── fig7_robustness.png
```

## Methodology

1. **Data preparation:** Use road-network node and road records derived from OpenStreetMap.
2. **Graph construction:** Represent road intersections or network locations as nodes and road connections as edges.
3. **Bridge analysis:** Identify bridge edges whose removal increases the number of connected components, using graph algorithms.
4. **Edge betweenness:** Analyze how frequently edges occur on shortest paths to identify structurally important connections.
5. **Link-addition scenarios:** Evaluate candidate additional links under defined budget constraints.
6. **Robustness analysis:** Compare network-level outcomes across the implemented scenarios.
7. **Visualization:** Generate maps and plots to communicate the findings.

## Figures and Results

The repository includes visualizations covering:

- Bridge locations in the road network
- Edge betweenness and structural impact
- EDI analysis
- Link-addition scenarios
- Candidate link augmentation
- Budget-related outcomes
- Network robustness

Open the PNG files in the repository to inspect the corresponding figures.

## How to Run

### 1. Clone the repository

```bash
git clone https://github.com/Navaraththinamlthursha/colombo-flood-resilience-graph-theory.git
cd colombo-flood-resilience-graph-theory
```

### 2. Install dependencies

```bash
pip install pandas numpy networkx matplotlib scipy
```

### 3. Run the analysis

```bash
python colombo_analysis.py
```

Execution time depends on the computer and the analysis performed. Check the script for input-file paths and output locations before running it.

The original exploratory notebook is also available as `notebook_original_colab.ipynb`.

## Data and Interpretation Limitations

- The analysis uses road-network structure rather than validated flood-event or flood-damage observations.
- Structural vulnerability should not be interpreted as a measured probability of flooding.
- Candidate link additions are analytical scenarios, not approved civil-engineering recommendations.
- Numerical results should be interpreted alongside the dataset definitions, preprocessing decisions, and assumptions in the analysis.

## Future Improvements

- Integrate authoritative flood-hazard and historical flood-event datasets.
- Evaluate accessibility to hospitals, emergency services, and evacuation destinations.
- Compare additional graph resilience metrics and alternative link-selection strategies.
- Validate proposed scenarios with transport and civil-engineering expertise.
- Improve reproducibility through documented data sources, dependency versions, and automated tests.

## Author

**Ithursha Navaraththinam**  
BSc (Hons) Data Science Undergraduate  
Sabaragamuwa University of Sri Lanka

## Acknowledgements

Road-network data source: OpenStreetMap. Please consult the applicable OpenStreetMap attribution and ODbL requirements when redistributing the data.

---

*This project was developed for academic research and learning in graph theory, network analysis, and infrastructure resilience.*

