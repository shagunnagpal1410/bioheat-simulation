# 🔥 Thermal Simulation with Finite Pointset Method (FPM)

This C++ project models **3D heat distribution** in biological tissue using the **Finite Pointset Method (FPM)**. It simulates temperature evolution over time due to three Gaussian needle heat sources. The model includes tissue and blood thermal properties and outputs the final temperature field.

---

## ⚙️ Running the Code

To run the simulation:

```
Ctrl + Shift + B
```

> Your `tasks.json` must be properly configured to compile and execute `6june.cpp`.  
> Ensure it includes the correct include path for the Eigen library.

---

## 📁 Required Input Files

Place these files in the same folder as `6june.cpp`:

| File Name             | Description                            |
|----------------------|----------------------------------------|
| `volume2_points.csv` | Volume mesh points (excluding tumor)   |
| `volume3_points.csv` | Tumor points to be merged with volume  |
| `boundary_nodes.csv` | Points on the tissue boundary (Dirichlet T=37°C) |

---

## 📦 Output Files

| File Name           | Description                                                    |
|--------------------|----------------------------------------------------------------|
| `6june.csv`        | Final coordinates and temperature values                        |
| `check_validity.csv` | Diagnostic info on unknown points and their interpolation weights |

---

## 🧠 Model Features

- Solves the **Pennes Bioheat Transfer Equation**
- Uses **Gaussian source model** for multiple needles
- Implements **meshfree Laplacian** via Moving Least Squares (MLS)
- Applies **Dirichlet boundary conditions** at the tissue surface
- Turns off heat sources when tissue temperature exceeds **105°C**
- Handles large 3D point clouds efficiently via voxelization

---

## 🧪 Applications

This simulation can be used for:

- **Thermal ablation modeling** in liver or tumor tissue
- Studying **bioheat propagation** in heterogeneous media
- Educational demonstrations of **meshfree methods (FPM)**

---

## 📁 Dependencies

- [Eigen Library](https://eigen.tuxfamily.org/)  
  Required for sparse matrix and linear algebra operations.

> Ensure Eigen is included in your include paths inside `tasks.json`.
