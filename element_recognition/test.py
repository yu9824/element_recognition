from element_recognition import get_ratio, make_compositions

if __name__ == '__main__':
    products = ['Li2La2O4']
    materials = ['Li2O', 'La2O3', 'TiO2']
    df = get_ratio(products, materials)
    print(df)

    print(make_compositions(materials))