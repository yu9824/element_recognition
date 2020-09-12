import numpy as np
import pandas as pd
import re

# from pdb import set_trace
from itertools import combinations
from copy import copy

default_elements = ("H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og")
# --- tools ---
def flatten_and_chomp(x):
    # とにかく一次元のリスト (or np.ndarray) 化
    if type(x) == np.ndarray:
        y = x.flatten()
    elif type(x) == str:
        y = [x]
    else:
        y = np.array(x).flatten()
    # 空白を削除してリストとして統一
    return list(map(lambda x:x.replace(' ',''), y))

'''
++ INPUT ++
compositions = 組成．一次元のリストを前提．(文字列も可)
elements: default "Og" (オガネソン, Oganesson)までの118元素のリスト．元素認識を行うにあたって必要な元素リストをカスタマイズ可能．

++ OUTPUT ++
pd.DataFrame, columns = elements, index = compositions, compositions中に含まれる元素の数をカウント．
'''
# --- main function ---
def ElementRecognition(compositions, **options):
    lst_compositions = flatten_and_chomp(compositions)

    elements = options['elements'] if 'elements' in options else default_elements
    delimiter = ('(', ')')  #'()'認識用

    # 単原子イオンと多原子イオンを判別し，別々のリストに分ける．
    monatomic_ion = []
    polyatomic_ion = []
    for ele in elements:
        if len(re.findall('[A-Z]', ele)) == 1:
            monatomic_ion.append(ele)
        else:
            polyatomic_ion.append(ele)


    ratio_ = np.zeros([len(lst_compositions), len(elements)])
    for i, comp in enumerate(lst_compositions):
        comp = lst_compositions[i]
        del_index = [j for j, x in enumerate(comp) if x.isupper() or x in delimiter]    # 大文字や()のときで区切り文字．

        dict_poly = {poly: comp.find(poly) for poly in polyatomic_ion if comp.find(poly) != -1}
        for poly, k in dict_poly.items():   # 多原子イオンを持っている場合は多原子イオンの部分を分割する値を削除．
            del_index = [j for j in del_index if j <= k or k + len(poly) <= j]

        dict_comp = {}  #元素の位置(index)と組成数(value)を取っておく辞書
        dict_brackets = {'f':[], 'b':[], 'value':[]}    # かっこの開始位置 (f)， 終了位置 (b), かっこ内の元素を何倍するか (value) を取っておく辞書
        for j in range(len(del_index)): #'(', ')', 'polyatomic_ion', 'monatomic_ion'がそれぞれどこか記録．
            # aとして抜き出す部分を取り出す．
            a = comp[del_index[j]:] if j == len(del_index) - 1 else comp[del_index[j]:del_index[j+1]]

            # 取り出したaを数値なのか，()なのかなど認識し，数字と元素に分けるなどする．
            if a in polyatomic_ion: #'NH4'(polyatomic_ion)
                dict_comp[a] = {'index':del_index[j], 'value':1.0}
            elif a == delimiter[0]: #'('
                dict_brackets['f'].append(del_index[j])
            elif a == delimiter[1]: #')'
                dict_brackets['b'].append(del_index[j])
                dict_brackets['value'].append(1.0)
            elif delimiter[1] in a: #')n'
                dict_brackets['b'].append(del_index[j])
                dict_brackets['value'].append(float(a.strip(delimiter[1])))
            else:   #monatomic_ion
                for k, s in enumerate(a):
                    if not s.isalpha():
                        dict_comp[a[:k]] = {'index':del_index[j], 'value':float(a[k:])}
                        break
                    elif k + 1 == len(a):   #最後までs.isalpha() == Trueとならなかった場合．
                        dict_comp[a] = {'index':del_index[j], 'value':1.0}

        #()のペアを一致させるための処理
        dict_brackets['f'].reverse()
        for j in range(len(dict_brackets['f'])):    #()のペアを一致させる
            if dict_brackets['f'][j] > dict_brackets['b'][j]:
                k = 1
                while j + k < len(dict_brackets['f']):
                    if dict_brackets['f'][j + k] < dict_brackets['b'][j]:
                        dict_brackets['f'][j], dict_brackets['f'][j + k] = dict_brackets['f'][j + k], dict_brackets['f'][j]
                        break
                    else:
                        k += 1

        for ele in dict_comp:
            for j in range(len(dict_brackets['f'])):    #()の中に入ってる元素の組成数をn倍にする
                if dict_brackets['f'][j] < dict_comp[ele]['index'] and dict_comp[ele]['index'] < dict_brackets['b'][j]:
                    dict_comp[ele]['value'] *= dict_brackets['value'][j]
            for j in range(len(elements)):  # 含まれてる元素のそのものの値を抜き出す．
                if ele == elements[j]:
                    ratio_[i][j] = dict_comp[ele]['value']

    df_output = pd.DataFrame(ratio_.copy(), index = lst_compositions, columns = elements)

    return df_output

'''
++ INPUT ++
products: 生成物，文字列，リストどちらでも．(必須)
materials: 原料．リスト．無駄な原料を入れると計算できなくなることが多いので入れるべきではない．(必須)
mathch_all: Default: False, boolean. 検算してすべての元素の割合が合ってるかを確かめ，一つでも元素の数が合わないと

++ OUTPUT ++
pd.DataFrame, columnsはmaterials, indexはproducts, それぞれの割合が入ってる．
もし負の割合がある場合は，その原料を混ぜただけではその生成物ができないことを表す．
'''

def Ratio(products, materials, **options):
    products = flatten_and_chomp(products)
    materials = flatten_and_chomp(materials)

    df_products = ElementRecognition(products)
    df_materials = ElementRecognition(materials)

    products_nonzero = df_products.iloc[:, list(set(df_products.values.nonzero()[1]))] # pd.Series
    materials_nonzero = df_materials.loc[:, products_nonzero.columns].transpose() # pd.Seriesに合わせてindexに元素 ('Li'etc.) が来るようにtranspose()

    counta_nonzero = np.count_nonzero(materials_nonzero, axis = 1)
    index_list = [list(np.where(counta_nonzero == n)[0]) for n in range(len(counta_nonzero), 0, -1) if len(np.where(counta_nonzero == n)[0]) != 0]  # 降順でmaterialsに共通で入ってる数が多い順に並べ， 同じ数共通で入っている場合は二次元ベクトルで表現
    i = 0

    # Ax = bのAを正方行列かつrank(A)=len(A)に整える．
    del_index_cand = [np.where(counta_nonzero == 0)[0]] # すべてが0 (materialsには入ってないけどproductsには入ってる元素) の組成比を削除． これはもし空集合であったとしても削除をしない，という条件を探れるので良い．　これを入れた理由は，基本的にnon_zeroの数が多い方から削除候補として用いているため，non_zero = 0が最後に削除されることになってしまうが， この条件は本来計算に含まれて胃はいけない事項であるから．
    indexes_memo = []
    for indexes in index_list:
        for i in range(1, len(indexes)):
            del_index_cand.extend(indexes_memo + list(map(list, combinations(indexes, i))))
        else:
            indexes_memo += indexes
            del_index_cand.append(copy(indexes_memo))

    if np.linalg.matrix_rank(materials_nonzero.to_numpy()) >= len(materials):
        for del_index in del_index_cand:
            A = np.delete(materials_nonzero.to_numpy(), del_index, axis = 0)
            if np.linalg.matrix_rank(A) == len(materials) and A.shape[0] == A.shape[1]:   # 正方行列で， またはそれと同じ大きさのrankを持っているとき
                break
        else:
            print(materials_nonzero, products_nonzero)
            exit("We can't solve.\nWe can't get square matrix.")
    else:
        print(materials_nonzero, products_nonzero)
        exit("We can't solve.\nThe rank(A) is lower than a number of variables(materials).")

    match_all = options['match_all'] if 'match_all' in options else False
    df_output = pd.DataFrame()
    for name, series in products_nonzero.iterrows():
        b = np.delete(series.to_numpy(), del_index)
        x = np.linalg.solve(A, b)
        sr_x = pd.Series(x, index = materials_nonzero.columns, name = name)
        if match_all:
            ar_memo = np.zeros(df_products.shape[1])
            for y, z in zip(x, df_materials.to_numpy()):
                ar_memo += y * z
            if np.allclose(ar_memo, df_products.loc[name].to_numpy()):    # すべての元素が検算で正しいとされるならば
                df_output = pd.concat([df_output, sr_x], axis = 1, sort = False)
        else:
            df_output = pd.concat([df_output, sr_x], axis = 1, sort = False)
    return df_output.transpose()




if __name__ == '__main__':
    pd.set_option('display.max_columns', 150)

    # products = 'Li0.33La0.5TiO3'
    products = ['Li2LaTiO6', 'Li0.33La0.5TiO3', 'Li2LaTiO6']
    materials = ['Li2O', 'LaO3', 'TiO2']

    # df_er = ElementRecognition(products)
    df_r = Ratio(products, materials, match_all = True)
    print(df_r)
    # print(df_er, df_r)
