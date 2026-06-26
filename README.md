# SCTOP

MATLAB implementation for paper:  
**"An Efficient Sub-Cell Topology Optimization Method: MATLAB Implementation for 3D Lattices"**  
*Chao Li, Qihan Wang, Minghui Zhang, Wei Gao, Zhen Luo (2026)*

## Usage

```matlab
SCTOP(nelx,nely,nelz,edgeN,volfrac,penal,rmin)
```

| Input | Description |
|-------|-------------|
| `nelx, nely, nelz` | Elements in x,y,z directions |
| `edgeN` | Sub-cells number along each direction of an element |
| `volfrac` | Maximum allowed volume fraction |
| `penal` | Penalization factor (SIMP, typically 3) |
| `rmin` | Filter radius |

## Example

```matlab
SCTOP(10,10,10,3,0.2,3,3);
```

---

**Last Update:** 26.06.2026