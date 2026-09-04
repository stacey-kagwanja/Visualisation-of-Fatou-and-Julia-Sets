# Visualisation-of-Fatou-and-Julia-Sets
*Work in progress*

This project explores the computer-graphical generation of Julia sets corresponding to various functions. 

---

## Overview

After continuously iterating a set of points through a rational function, subsets of these points called *basins* converge to either infinity or a point on the complex plane. The shared (topological) boundary of these points is called the Julia set, and the union of the basins form the Fatou set corresponding to the function.

The points in each basin can be coloured according to how many iterations it takes for points to approach sufficiently near whichever point they converge to.

Further details of this are in my essay *Fatou and Julia Sets of Holomorphic Functions on the Riemann Sphere*

---

## Project Components

---

## Tech Stack

| Component | Technology |
|-----------|------------|
| Language | Python |

---

## Installation

Clone the repository and install the required dependencies:
```bash
git clone https://github.com/stacey-kagwanja/Visualisation-of-Fatou-and-Julia-Sets.git
cd Visualisation-of-Fatou-and-Julia-Sets
pip install -r requirements.txt
```

---

## Usage

Download the Jupyter notebook (`Fatou and Julia.py`), then open and run the notebook.

```
Fatou and Julia.py
```

Running this code generates an example image of the Fatou and Julia sets, of the function z^2 + 1.


Download the Jupyter notebook (`Fatou and Julia.py`), then open and run the notebook.

```
Fatou and Julia Images.ipynb
```

Running this code generates example images of the Fatou and Julia sets of some rational functions. It also includes code for creating colour maps to use for each basin.

---

## Project Context

...

---

## Author

- Stacey Kagwanja

---

## Acknowledgements

- ...
