# element_recognition
文字列として書かれた組成をDataFrameにして整理したり，逆に混合比率から組成を生成するパッケージ

## 使い方
### モジュールインポート
```python
from element_recognition import get_ratio, make_compositions
```

#### ```get_ratio```
```python
get_ratio(products = ['LiLaTiO6'], materials = ['Li2O', 'LaO3', 'TiO2'])
```
