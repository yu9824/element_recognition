# Element Recognition
![python_badge](https://img.shields.io/pypi/pyversions/element-recognition)
![license_badge](https://img.shields.io/pypi/l/element-recognition)
![Total_Downloads_badge](https://pepy.tech/badge/element-recognition)

## What is this.
A package that organizes compositions written as strings into pandas.DataFrames, or conversely, generates compositions from mixing ratios by using only `numpy` and `pandas`.

It works well with [ternary_diagram](https://pypi.org/project/ternary-diagram/) (This makes it easy to generate beautiful ternary diagram.).
Please check it out.

## How to install 
```bash
pip install element-recognition
```

## How to use
If you want to know in details, see [example](https://github.com/yu-9824/element_recognition/tree/main/example).


### Import modules
```python
from element_recognition import get_ratio, make_compositions
```

### ```get_ratio```
This is a function that returns a pandas.DataFrame of mixing ratios given a compound and a raw material.
```python
get_ratio(products = ['LiLa2TiO6'], materials = ['Li2O', 'La2O3', 'TiO2'])
```
```
               Li2O  La2O3  TiO2
    Li2La2TiO6   1.0   1.0   1.0
```

### ```make_compositions```
This function returns the composition formula and the amount of all elements contained as a pandas.DataFrame by giving the raw materials and the mixing ratio.
```python
make_compositions(materials = ['Li2O', 'La2O3', 'TiO2'], ratio = [[1, 2, 3]])
```
```
                H   He   Li   Be    B    C    N     O    F   Ne   Na   Mg   Al   Si    P  ...   Rf   Db   Sg   Bh   Hs   Mt   Ds   Rg   Cn   Nh   Fl   Mc   Lv   Ts   Og
    Li2Ti3La4O13  0.0  0.0  2.0  0.0  0.0  0.0  0.0  13.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  ...  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```

## LICENSE
see [LICENSE](https://github.com/yu-9824/element_recognition/tree/main/LICENSE).