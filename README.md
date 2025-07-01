# 🔥 Thermal Simulation with Finite Pointset Method (FPM)

This C++ project models **3D heat distribution** in biological tissue using the **Finite Pointset Method (FPM)**. It simulates temperature evolution over time due to three Gaussian needle heat sources. The model includes tissue and blood thermal properties and outputs the final temperature field.

---

## 📁 Input Files Required

Make sure the following CSV files are in the **same folder** as `6june.cpp` before running:

| File Name             | Description                            |
|----------------------|----------------------------------------|
| `volume2_points.csv` | Volume mesh points (excluding tumor)   |
| `volume3_points.csv` | Tumor points to be merged with volume  |
| `boundary_nodes.csv` | Points on the tissue boundary (Dirichlet T=37°C) |

---

## ⚙️ Dependencies

- **C++17 or later**
- [Eigen Library](https://eigen.tuxfamily.org/) (used for sparse linear algebra)

Make sure the Eigen headers are accessible during compilation.

---

## 🚀 How to Compile and Run

### Compile
```bash
g++ 6june.cpp -o thermal_sim -O2 -std=c++17
```

> Replace `g++` with `g++-11` or similar if required.

### Run
```bash
./thermal_sim
```

You will be prompted to enter:
- Maximum simulation time
- Time step size
- Coordinates and `alpha` values for three heat-delivering needles

---

## 📦 Output Files

| File Name         | Description                                             |
|------------------|---------------------------------------------------------|
| `6june.csv`       | Contains final coordinates and corresponding temperatures |
| `check_validity.csv` | Debug info: identifies known/unknown nodes and interpolation weights |

---

## 🧠 Model Features

- Solves the **Pennes Bioheat Transfer Equation**
- Uses **Gaussian source model** for multiple needles
- Implements **meshfree Laplacian** via Moving Least Squares (MLS)
- Applies **Dirichlet boundary conditions** at surface
- Disables heat sources after tissue exceeds **105°C**
- Efficiently handles large 3D point clouds with voxel-based spatial search

---

## 🧪 Applications

This simulation can be used for:

- **Thermal ablation modeling** in liver or tumor tissue
- Studying **bioheat propagation** in heterogeneous media
- Educational demonstrations of **meshfree methods (FPM)**
