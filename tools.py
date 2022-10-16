def indexing_values(liste, value):
    index = []
    for i in range(len(liste)):
        if type(liste[i]) == list:
            for j in range(len(liste[i])):
                if liste[i][j] == value:
                    index.append([i, j])
                elif type(liste[i][j]) == list:
                    for m in range(len(liste[i][j])):
                        if liste[i][j][m] == value:
                            index.append([i, j, m])
                        elif type(liste[i][j][m]) == list:
                            for l in range(len(liste[i][j][m])):
                                if liste[i][j][m][l] == value:
                                    index.append([i, j, m, l])
                                elif type(liste[i][j][m][l]) == list:
                                    for n in range(len(liste[i][j][m][l])):
                                        if liste[i][j][m][l][n] == value:
                                            index.append([i, j, m, l, n])
                                        elif type(liste[i][j][m][l][n]) == list:
                                            for o in range(len(liste[i][j][m][l][n])):
                                                if liste[i][j][m][l][n][o] == value:
                                                    index.append([i, j, m, l, n, o])
        elif type(liste[i]) == str:
            for k in range(len(liste[i])):
                if liste[i][k] == value:
                    index.append([i, k])
                elif type(liste[i][k]) == list:
                    for l in range(len(liste[i][k])):
                        if liste[i][k][l] == value:
                            index.append([i, k, l])
                        elif type(liste[i][k][l]) == list:
                            for m in range(len(liste[i][k][l])):
                                if liste[i][k][l][m] == value:
                                    index.append([i, k, l, m])
                                elif type(liste[i][k][l][m]) == list:
                                    for n in range(len(liste[i][k][l][m])):
                                        if liste[i][k][l][m][n] == value:
                                            index.append([i, j, m, l, n])
                                        elif type(liste[i][j][m][l][n]) == list:
                                            for o in range(len(liste[i][k][l][m][n])):
                                                if liste[i][k][l][m][n][o] == value:
                                                    index.append([i, j, m, l, n, o])

    return index

def get_element(index, tab=list):
    for ind in range(len(index)):
        try:
            i = index[ind][0]
        except: i = -1

        try:
            j = index[ind][1]
        except: j = -1

        try:
            k = index[ind][2]
        except: k = -1

        try:
            l = index[ind][3]
        except: l = -1

        try:
            m = index[ind][4]
        except: m = -1

        try:
            n = index[ind][5]
        except: n = -1

        try:
            o = index[ind][6]
        except: o = -1

        if i != -1:
            if j != -1:
                if k != -1:
                    if l != -1:
                        if m != -1:
                            if n != -1:
                                if o != -1: return tab[i][j][k][l][m][n][o]
                            else: return tab[i][j][k][l][m][n]
                        else: return tab[i][j][k][l][m]
                    else: return tab[i][j][k][l]
                else: return tab[i][j][k]
            else: return tab[i][j]
        else: return tab[i]

def diviser_liste(liste, divided: int):
    nouvelle_liste = []
    for i in range(0, len(liste), divided):
        nouvelle_liste.append(liste[i:i+divided])
    return nouvelle_liste

def diviser_elements(liste, number: int):
    nouvelle_liste = []
    for element in range(number):
        for sous_liste in liste:
            try:
                nouvelle_liste.append(sous_liste[element])
            except: pass
    return nouvelle_liste

