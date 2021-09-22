'''
Copyright © 2021 yu9824
'''

from collections import Counter
import numpy as np
import pandas as pd
import re

from itertools import combinations, product
from copy import copy
from math import ceil, floor


DEFAULT_ELEMENTS = ("H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og")
# --- tools ---
def _check_value(x:str) -> float:
    """check value of the elements

    Parameters
    ----------
    x : str

    Examples
    ----------
    >>> _check_value('1.2')
    1.2
    >>> _check_value('1/2')
    0.5

    Returns
    -------
    float
        checked value
    """
    c = Counter(x)
    if '/' in c and c['/'] == 1:
        _a, _b = map(float, x.split('/'))
        return _a / _b
    else:
        return float(x)

def _flatten_and_trim(x):
    # とにかく一次元のリスト (or np.ndarray) 化
    if isinstance(x, str):
        x = [x]
    _x = np.array(x).ravel()
    # 空白を削除してリストとして統一
    return list(map(lambda x:x.replace(' ',''), _x))

# --- main function ---
def element_recognition(compositions, elements = None) -> pd.DataFrame:
    '''
    Parameters
    ----------
    compositions: str or list
        組成．一次元のリストを前提．(文字列も可)
    elements: list, default None. It means "Og" (オガネソン, Oganesson)までの118元素のリスト．
        元素認識を行うにあたって必要な元素リストをカスタマイズ可能．

    Returns
    ----------
    pd.DataFrame, columns = elements, index = compositions, compositions中に含まれる元素の数をカウント．
    '''

    if elements is None:
        elements = DEFAULT_ELEMENTS

    lst_compositions = _flatten_and_trim(compositions)

    delimiter = ('(', ')')  #'()'認識用

    # 単原子イオンと多原子イオンを判別し，別々のリストに分ける．
    monatomic_ion = []
    polyatomic_ion = []
    for ele in elements:
        if len(re.findall('[A-Z]', ele)) == 1:
            monatomic_ion.append(ele)
        else:
            polyatomic_ion.append(ele)


    ratio_ = np.zeros([len(lst_compositions), len(elements)], dtype=float)
    for i, comp in enumerate(lst_compositions):
        del_index = [j for j, x in enumerate(comp) if x.isupper() or x in delimiter]    # 大文字や()のときで区切り文字．

        dict_poly = {poly: comp.find(poly) for poly in polyatomic_ion if comp.find(poly) != -1}
        for poly, k in dict_poly.items():   # 多原子イオンを持っている場合は多原子イオンの部分を分割する値を削除．
            del_index = [j for j in del_index if j <= k or k + len(poly) <= j]

        dict_comp = {}  # key: index(元素の位置), value: value(組成数), symbol(元素記号, 多原子イオンも含む)
        dict_brackets = {
            'f':[],     # かっこの開始位置
            'b':[],     # かっこの終了位置
            'value':[]  # かっこ内の元素を何倍するか
        }
        for j in range(len(del_index)): #'(', ')', 'polyatomic_ion', 'monatomic_ion'がそれぞれどこか記録．
            # aとして抜き出す部分を取り出す．
            a = comp[del_index[j]:] if j == len(del_index) - 1 else comp[del_index[j]:del_index[j+1]]

            # 取り出したaを数値なのか，()なのかなど認識し，数字と元素に分けるなどする．
            if a in polyatomic_ion: #'NH4'(polyatomic_ion)
                dict_comp[del_index[j]] = {
                    'value':1.0,
                    'symbol':a,
                }
            elif a == delimiter[0]: #'('
                dict_brackets['f'].append(del_index[j])
            elif a == delimiter[1]: #')'
                dict_brackets['b'].append(del_index[j])
                dict_brackets['value'].append(1.0)
            elif delimiter[1] in a: #')n'
                dict_brackets['b'].append(del_index[j])
                dict_brackets['value'].append(_check_value(a.strip(delimiter[1])))
            else:   #monatomic_ion
                for k, s in enumerate(a):
                    if not s.isalpha():
                        dict_comp[del_index[j]] = {
                            'value':_check_value(a[k:]),
                            'symbol':a[:k],
                        }
                        break
                else:   #最後までs.isalpha() == Falseとならなかった (数字が現れなかった) 場合．
                    dict_comp[del_index[j]] = {
                        'value':1.0,
                        'symbol':a,
                    }

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
        
        for index in dict_comp:
            for j in range(len(dict_brackets['f'])):    #()の中に入ってる元素の組成数をn倍にする
                if dict_brackets['f'][j] < index and index < dict_brackets['b'][j]:
                    dict_comp[index]['value'] *= dict_brackets['value'][j]
            for j in range(len(elements)):  # 含まれてる元素のそのものの値を抜き出す．
                if dict_comp[index]['symbol'] == elements[j]:
                    ratio_[i][j] += dict_comp[index]['value']

    df_output = pd.DataFrame(ratio_.copy(), index = lst_compositions, columns = elements)

    return df_output


def get_ratio(products, materials, exact = True, elements = None) -> pd.DataFrame:
    '''
    Parameters
    ----------
    products: str or list
        生成物，文字列，リストどちらでも．(必須)
    
    materials: list
        原料．リスト．無駄な原料を入れると計算できなくなることが多いので入れるべきではない．(必須)
    
    exact: bool, default True
        検算してすべての元素の割合が合ってるかを確かめ，一つでも元素の数が合わないとき，Noneとして返す．

    elements: list, default None.
        This is variable for element_recognition function.

    Returns
    ----------
    pd.DataFrame, columnsはmaterials, indexはproducts, それぞれの割合が入ってる．
    もし負の割合がある場合は，その原料を混ぜただけではその生成物ができないことを表す．

    Examples
    ----------
    >>> materials = ['Li2O', 'La2O3', 'TiO2']
    >>> products = ['Li2La2TiO6']
    >>> get_ratio(products, materials)
                Li2O  La2O3  TiO2
    Li2La2TiO6   1.0    1.0   1.0
    '''
    products = _flatten_and_trim(products)
    materials = _flatten_and_trim(materials)

    df_products = element_recognition(products)
    df_materials = element_recognition(materials)

    products_nonzero = df_products.iloc[:, list(set(df_products.values.nonzero()[1]))] # pd.Series
    materials_nonzero = df_materials.loc[:, products_nonzero.columns].transpose() # pd.Seriesに合わせてindexに元素 ('Li'etc.) が来るようにtranspose()

    counta_nonzero = np.count_nonzero(materials_nonzero, axis = 1)
    index_list = [list(np.where(counta_nonzero == n)[0]) for n in range(len(counta_nonzero), 0, -1) if len(np.where(counta_nonzero == n)[0]) != 0]  # 降順でmaterialsに共通で入ってる数が多い順に並べ， 同じ数共通で入っている場合は二次元ベクトルで表現

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
            print(materials_nonzero, '\n', products_nonzero)
            raise ValueError("We can't solve.\nWe can't get square matrix.")
    else:
        raise ValueError("We can't solve.\nThe rank(A) is lower than a number of variables(materials).")

    df_output = pd.DataFrame()
    _srs = []
    for name, series in products_nonzero.iterrows():
        b = np.delete(series.to_numpy(), del_index)
        try:
            x = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            raise np.linalg.LinAlgError("We can't solve.\n", A, b)
        sr_x = pd.Series(x, index = materials_nonzero.columns, name = name)
        if exact:
            ar_memo = np.zeros(df_products.shape[1])
            for y, z in zip(x, df_materials.to_numpy()):
                ar_memo += y * z
            if np.allclose(ar_memo, df_products.loc[name].to_numpy()):    # すべての元素が検算で正しいとされるならば
                # df_output = pd.concat([df_output, sr_x], axis = 1, sort = False)
                _srs.append(sr_x)
            else:
                _srs.append(pd.Series([None] * len(materials_nonzero.columns), index = materials_nonzero.columns, name = name))
                # df_output = pd.concat([df_output, pd.Series([None] * len(materials_nonzero.columns), index = materials_nonzero.columns, name = name)], axis = 1, sort = False)
        else:
            _srs.append(sr_x)
            # df_output = pd.concat([df_output, sr_x], axis = 1, sort = False)
    df_output = pd.concat(_srs, axis=1).transpose()
    return df_output


def make_compositions(materials, ratio = None, easy = True, max_comp:int = 15, front = None, back = None, max_show_prec:int = 3, elements = None) -> pd.DataFrame:
    '''
    Parameters
    ----------
    materials : list
        e.g.) ['Li2O', 'La2O3', 'TiO2']

    ratio: 2D-list (pd.DataFrame, np.ndarray etc.), default None
        If no input is given, an appropriate composition is generated.

    easy: bool, default True
        Whether or not to lighten the density of the composition generated when the ratio is None; lighten when True.

    max_comp: int, default 15
        The maximum number of composition ratios that will be automatically generated  when the ratio is None.

    front: list, default None. It means ('Li', 'Na', 'K', 'Rb', 'Cs'). 
        An element that is preferentially in front of a composition when it is generated. The more elements are in front, the more priority is given to the front. If duplicate elements are specified, back is given priority.
    
    back: list, default None. It means ('I', 'Br', 'Cl', 'F', 'S', 'O').
        An element that is preferentially placed behind a composition when it is generated. The more elements are in front, the more elements are in the back. If duplicate elements are specified, back is given priority.

    max_show_prec: int, default 3
        The maximum number of float places to display when creating a compound.

    elements: list, default None.
        This is variable for element_recognition function.

    Returns
    ----------
    pd.DataFrame. The index represents the composition of the generated material.
    '''

    # frontとbackに関する処理
    default_front = ('Li', 'Na', 'K', 'Rb', 'Cs')
    default_back = ('I', 'Br', 'Cl', 'F', 'S', 'O')
    # 重複を削除
    if front is None and back is None:
        front = default_front
        back = default_back
    elif front is None:
        front = [ele for ele in default_front if ele not in back]
    elif back is None:
        back = [ele for ele in default_back if ele not in front]
    else:   # 両方が指定された場合はbackを優先
        front = [ele for ele in front if ele not in back]

    # 何も入力されなければ自動生成
    p = product(np.arange(max_comp+1), repeat = len(materials))
    if ratio is None:
        ratio = [i for i in p if sum(i) == max_comp] if easy else [i for i in p if sum(i) != 0]
    ratio = np.array(ratio)

    if len(ratio.shape) == 1 and ratio.shape[0] == len(materials):
        ratio = ratio.reshape(1, -1)
    elif len(ratio.shape) != 2 or ratio.shape[1] != len(materials):
        raise ValueError('A shape of ratio is not correct; The shape is', ratio.shape)

    df_materials = element_recognition(materials)
    df_products = pd.DataFrame(ratio @ df_materials.to_numpy(), columns = df_materials.columns)

    compositions = []
    for row in df_products.to_numpy():
        dict_memo = {}  # keyが元素名，valueがたすべき文字列の辞書を作成
        for i, x in enumerate(row):
            if x == 0:
                pass
            elif x == 1:
                dict_memo[df_materials.columns[i]] = df_materials.columns[i]
            elif ceil(x) == floor(x):   # 0と1以外の整数のとき
                dict_memo[df_materials.columns[i]] = df_materials.columns[i] + str(ceil(x))
            else:   # 小数のとき
                _temp = ['{0:.{1}f}'.format(x, prec) for prec in range(max_show_prec+1)]
                base = _temp.pop(-1)
                s = base
                while _temp:
                    cand = _temp.pop(-1)
                    if float(cand) == float(base):
                        s = cand
                    else:
                        break
                dict_memo[df_materials.columns[i]] = df_materials.columns[i] + str(s)

        # 本当はfront, backが正しく入力されてるか考えなきゃいけなさそう．
        # 最後joinして文字列にするmemoリストを定義
        memo = []

        # frontに入ってるのを優先的に取り出す． (追加したら削除)
        for f in front:
            if f in dict_memo:
                memo.append(dict_memo.pop(f))
        # backに入ってないのを追加していく． (追加したら削除)
        keys = dict_memo.copy().keys()
        for k in keys:
            if k not in back:
                memo.append(dict_memo.pop(k))
        # 残ったbackの逆順に追加．(追加したら削除)
        for b in reversed(back):
            if b in dict_memo:
                memo.append(dict_memo.pop(b))
                if len(dict_memo) == 0:
                    break

        compositions.append(''.join(memo))

    df_products.index = compositions
    return df_products



if __name__ == '__main__':
    # from pdb import set_trace
    # pd.set_option('display.max_columns', 150)

    # products = ['LaO3', 'Li2LaO']
    # # products = 'Li0.33La0.5TiO3'
    # products = ['Li2LaTiO6', 'Li0.33La0.5TiO3', 'Li2LaTiO6']
    # # products = ['Li2LaTiO6', 'Li0.33La0.5TiO3', 'Li2LaTiO6', 'Li2O']
    # materials = ['Li2O', 'LaO3', 'TiO2']

    # # products = 'LiLaO2'
    # # products = ['Li2La2O4', 'LiLaO2']
    # # materials = ['Li2TiO3', 'Li2BaO3', 'LiLaO2']

    # df_er = element_recognition(products)
    # df_r = Ratio(products, materials, exact = True)
    # print(df_r)
    # # print(df_er, df_r)

    # df_m = make_compositions(materials, ratio = [2/3, 2/3, 3])
    # print(df_m)
    # # print(df_m)
    from doctest import testmod
    testmod(verbose=True)
