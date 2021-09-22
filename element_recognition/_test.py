if __name__ == '__main__':
    from element_recognition import get_ratio, make_compositions
    import numpy as np
    import pandas as pd

    products = ['Li2La2O4', 'Li0.5La1/2TiO3', '(LiLa)1/2TiO3', 'Li2O']
    materials = ['Li2O', 'La2O3', 'TiO2']
    df = get_ratio(products, materials)
    print(df)

    df = get_ratio(products[0], materials)
    print(df)

    df = get_ratio(np.array(products), materials)
    print(df)

    df = get_ratio(pd.Series(products), materials)
    print(df)

    df = get_ratio(pd.DataFrame(products), materials)
    print(df)

    print(make_compositions(materials))