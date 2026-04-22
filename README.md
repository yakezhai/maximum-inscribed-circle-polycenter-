
# maximum-inscribed-circle-polycenter: Fast and Precise Polygon Center Identification

**Polycenter** is a high-performance, exact algorithm for finding the maximum inscribed circle (visual center) of any non-self-intersecting polygon, including those with holes and parallel sides. 

Unlike approximation methods like *Polylabel*, Polycenter achieves **global optimality** and operates at the limit of floating-point precision without requiring a predefined precision parameter.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Paper: IJGIS](https://img.shields.io/badge/Journal-IJGIS-red.svg)](https://doi.org/10.1080/13658816.2025.2514056)

---

## Key Features

- **Mathematically Exact**: Combines binary search for local optima with a divide-and-conquer strategy for global optimality.
- **Superior Performance**: 
  - **3x faster** than Polylabel at a precision of $1.0$.
  - **10x faster** than Polylabel at a precision of $10^{-8}$.
- **Robustness**: Handles complex geometries, islands, holes, and parallel sides (where Polylabel often fails).
- **Header-Only**: Single `.hpp` file. No external dependencies. Easy to integrate into any C++ project.
- **Precision-Independent**: Leverages polygon elements in contact with the circle, eliminating the need for user-defined tolerance.

## Applications

Polycenter is designed for diverse fields requiring spatial analysis:
- **GIS**: Optimal label placement.
- **Robotics**: Obstacle avoidance and object grasping.
- **Materials Science**: Roundness and sphericity evaluation.
- **Urban Planning & Image Recognition**.

---

## Citation

If you use this algorithm or the provided data in your research, please cite our official **IJGIS** paper:

```bibtex
@article{ijgis2026_polycenter,
  author = {Yake Zhai},
  title = {Polycenter: fast and precise polygon center identification for GIS and beyond},
  journal = {International Journal of Geographical Information Science},
  volume = {40},
  number = {2},
  pages = {327--347},
  year = {2026},
  publisher = {Taylor & Francis},
  doi = {10.1080/13658816.2025.2514056},
  url = {[https://doi.org/10.1080/13658816.2025.2514056](https://doi.org/10.1080/13658816.2025.2514056)}
}

figshare:  https://figshare.com/articles/dataset/Code_and_test_data_for_Polycenter/28244642
