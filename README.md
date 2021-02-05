# Element Recognition
A package that organizes compositions written as strings into pandas.DataFrames, or conversely, generates compositions from mixing ratios by using only `numpy` and `pandas`.

It works well with [ternary_diagram](https://github.com/yu-9824/ternary_diagram) (a self-made package for easy phase diagram generation) if you want, so please check it out.

## How to use
詳細は[example](https://github.com/yu-9824/element_recognition/tree/master/example)を参照．

### Import modules
```python
from element_recognition import get_ratio, make_compositions
```

#### ```get_ratio```
```python
get_ratio(products = ['LiLaTiO6'], materials = ['Li2O', 'LaO3', 'TiO2'])
```
```
               Li2O  LaO3  TiO2
    Li2LaTiO6   1.0   1.0   1.0
```

#### ```make_compositions```
```python
make_compositions(materials = ['Li2O', 'LaO3', 'TiO2'], ratio = [[1, 2, 3]])
```
```
                H   He   Li   Be    B    C    N     O    F   Ne   Na   Mg   Al   Si    P  ...   Rf   Db   Sg   Bh   Hs   Mt   Ds   Rg   Cn   Nh   Fl   Mc   Lv   Ts   Og
    Li2Ti3La2O13  0.0  0.0  2.0  0.0  0.0  0.0  0.0  13.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  ...  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```