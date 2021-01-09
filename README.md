# element_recognition
文字列として書かれた組成をDataFrameにして整理したり，逆に混合比率から組成を生成するパッケージ

## 使い方
詳細はexampleを参照．
### モジュールインポート
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
make_cimpositions(materials = ['Li2O', 'LaO3', 'TiO2'], ratio = [[1, 2, 3]])
```
```
                H   He   Li   Be    B    C    N     O    F   Ne   Na   Mg   Al   Si    P  ...   Rf   Db   Sg   Bh   Hs   Mt   Ds   Rg   Cn   Nh   Fl   Mc   Lv   Ts   Og
    Li2Ti3La2O13  0.0  0.0  2.0  0.0  0.0  0.0  0.0  13.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  ...  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```