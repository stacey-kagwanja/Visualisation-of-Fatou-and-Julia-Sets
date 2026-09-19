# Visualisation-of-Fatou-and-Julia-Sets

This project explores the computer-graphical generation of Julia sets corresponding to various functions. 

---

## Overview

After continuously iterating a set of points through a rational function, subsets of these points called *basins* converge to either infinity or a point on the complex plane. The shared (topological) boundary of these points is called the Julia set, and the union of the basins form the Fatou set corresponding to the function.

The points in each basin can be coloured according to how many iterations it takes for points to approach sufficiently near whichever point they converge to.

Further details of this are in my essay *Fatou and Julia Sets of Holomorphic Functions on the Riemann Sphere*.

<img src="Images%20for%20README/Fatou%20and%20Julia%20sets%20for%20(2z3%20+%201)%25(3z2)%20image.png" width="600"/>

<img src="Images%20for%20README/Fatou%20and%20Julia%20sets%20for%20(2z3%20+%201)%25(3z2)%20zoom%20in%20GIF.gif" width="600"/>

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


Download the Jupyter notebook (`Fatou and Julia Images.ipynb`), then open and run the notebook.

```
Fatou and Julia Images.ipynb
```

Running this code generates example images of the Fatou and Julia sets of some rational functions. It also includes code for creating colour maps to use for each basin.

---

## Project Context

The code from this project was used to generate images and GIFs for my university essay *Fatou and Julia Sets of Holomorphic Functions on the Riemann Sphere* and the corresponding presentation.

---

## Author

- Stacey Kagwanja

---

## Acknowledgements

- The code was adapted from the below source to work for any holomorphic function on the Riemann sphere, to illustrate the basins of the attracting periodic orbits by colouring points in each basin based on the number of iterations taken for each point to arrive within some tolerance of the corresponding orbit, and to reduce the runtime for GIF generation by producing several frames at a time:
  
  B. L. Badger. Julia sets — form and formula: https://blbadger.github.io/julia-sets.html
