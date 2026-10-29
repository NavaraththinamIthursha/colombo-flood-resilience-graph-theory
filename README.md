# Colombo Road Network - Graph Theory Research Assignment

Title: Strengthening the Colombo Road Network: A Graph-Theoretic Study of Bridges,
Edge Betweenness and Budget-Constrained Link Addition

## Contents
| File / folder | What it is |
|---|---|
| `Colombo_Graph_Theory_Report` (.docx / .pdf) | The research report (export from the online document) |
| `colombo_analysis.py` | One script that reproduces every number, table and figure |
| `colombo_nodes.csv`, `colombo_roads.csv` | OpenStreetMap extract (1,737 nodes, 3,633 road records) |
| `figures/` | Figures 1-7 used in the report (PNG) |
| `results/` | Result tables (CSV) and `key_numbers.json` |
| `notebook_original_colab.ipynb` | Original Colab notebook (exploratory work) |
| `assignment_brief.pdf` | Assignment instructions |

## How to reproduce
    pip install pandas numpy networkx matplotlib scipy
    python colombo_analysis.py     # about 2-4 minutes; writes figures/ and results/

Random seeds are fixed (42, and 7 for the single random link set).

## Before you submit (checklist)
- [ ] Put your name, student ID and lecturer on the report title page
- [ ] Complete the 3 references marked [complete ... from the paper] (Kalpana 2021, Rebally 2021, He 2024, Boeing & Ha 2024 journal details) by opening each paper
- [ ] Read the four papers' full text and adjust Section 2 if anything differs
- [ ] Be ready to explain: bridge, Tarjan algorithm, edge betweenness, ILP set-cover formulation, why greedy leaves 8 bridges
- [ ] Meet your lecturer (the brief encourages this) and state the flood-data limitation
- [ ] Export report to PDF / Word as your department requires

## Honest notes
- No flood layer was available, so results are structural scenarios, not flood-damage estimates.
- 357 of 398 bridges are dead-end stubs; the link plan is a mathematical exercise, not an engineering proposal.
- Raw file gives 2,318 undirected pairs (one self-loop); the cleaned graph used for all results has 2,317 edges.
