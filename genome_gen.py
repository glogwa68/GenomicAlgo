import random
import numpy as np
import base64
import string
import threading
import tools
import tools as tool
import queue
from numba import cuda

global recurent
recurent = 0
celulle_count = 0
threadLimiter = threading.BoundedSemaphore(100)

def threaded(f, daemon=False):
    import queue

    q = queue.Queue(maxsize=100)

    def wrapped_f(q, *args, **kwargs):
        '''this function calls the decorated function and puts the
        result in a queue'''
        ret = f(*args, **kwargs)
        q.put(ret)

    def wrap(*args, **kwargs):
        '''this is the function returned from the decorator. It fires off
        wrapped_f in a new thread and returns the thread object with
        the result queue attached'''

        q = queue.Queue()

        t = threading.Thread(target=wrapped_f, args=(q,)+args, kwargs=kwargs)
        t.daemon = daemon
        t.start()
        t.result_queue = q
        return t

    return wrap

class cellule:
    def membrane():
        #gen taille fixe de membrane pour cellule
        return random.choice([1000, 2000, 4000, 8000, 16000, 32000, 64000])

    def noyau():
        #chromosome
        chromosome_final = chromosome.M_chromosome()
        return tool.diviser_liste(chromosome_final, 4)

    def M_cellule():
        #noyau
        #membrane
        chromosome_final = cellule.noyau()
        chromosome_final.append(cellule.membrane())

        return chromosome_final

class chromosome:
    def generate_double_chormatides():
        #chromatides x2 (assemblé)
        chromatides_simple = chromatides.M_chromatides()
        chromatides_double = tool.diviser_liste(chromatides_simple, 2)

        return chromatides_double

    def M_chromosome():
        #generate double chromatides add to tab 2 +100 a chaque telomere
        double_chromatides = chromosome.generate_double_chormatides()
        double_chromatides_scinder = tool.diviser_elements(double_chromatides, 4)
        chromosome_final = tool.diviser_liste(double_chromatides_scinder, random.choice([1000, 2000, 4000, 8000, 16000, 32000, 64000]))

        return chromosome_final

class chromatides:
    def telomere():
        #double brins ADN
        adn_complet = nucléosomes.M_nucléosomes()
        adn_scinder = tool.diviser_liste(adn_complet, 4)

        return adn_scinder

    def centromere():
        #connect fin de brin telomere to exact positon de debut telomere inversé superieur x4
        adn_scinder = chromatides.telomere()
        chromosome = tool.diviser_elements(adn_scinder, 4)

        return chromosome

    def M_chromatides():
        #generate 2 tab with 4 conection after 1er=> aux mileur a droite, 2eme=> milieux gauche connect 4 branche of 1/2 per tab
        chromosome = chromatides.centromere()
        chromatides_pass1 = tool.diviser_liste(chromosome, 2)
        chromatides_pass2 = tool.diviser_liste(chromatides_pass1, 2)

        chromatides_resort = []
        for chromosome in chromatides_pass2:
            chromatides_resort.append(chromosome)

        return chromatides_resort

class nucléosomes:
    def calcul_adn(division_int: int):
        HPON = adn.get_acide_HPON()

        H_sum = 0
        P_sum = 0
        O_sum = 0
        N_sum = 0
        for H_P_O_N in range(len(HPON)):
            for mol in range(len(HPON[H_P_O_N])):
                if H_P_O_N == 0:
                    H_sum += sum(HPON[H_P_O_N][mol])
                elif H_P_O_N == 1:
                    P_sum += sum(HPON[H_P_O_N][mol])
                elif H_P_O_N == 2:
                    O_sum += sum(HPON[H_P_O_N][mol])
                elif H_P_O_N == 3:
                    N_sum += sum(HPON[H_P_O_N][mol])

        adn_fullpart = adn.M_adn(int(division_int), H_sum, P_sum, O_sum, N_sum)  # range_of_adn
        part_adn = []
        for molecule in adn_fullpart:
            divided_molecule = tool.diviser_liste(molecule, random.randint(int(4*8198), int(4*65584)))
            parted = 0
            for part in divided_molecule:
                for i in range(len(part)):
                    parted += int(part[i])
            part_adn.append(parted)

        sum_adn = []
        divided_adn = tools.diviser_liste(part_adn, random.randint(10*2, 10*16))
        for adn_part in divided_adn:
            sum_adn.append(sum(adn_part))

        print(f"(cell divided by {division_int}) :: resulted list of {len(sum_adn)} cells.")
        return sum_adn

    def histones():
        #enroule une elise d'adn en bobine de 1 lise d'adn de longeur d'un chaine d'adn et les fusione avec chaque brin.
        #prend en entré l'adn pour le repliqué pour donné des nucléosomes
        range_of_adn = random.randint(100000000, 180000000)

        master_adn = nucléosomes.calcul_adn(range_of_adn)
        formed_adn = []
        for summing_adn in range(len(master_adn)):
            sum1 = random.choice(master_adn)
            sum2 = random.choice(master_adn)
            if sum1 != sum2:
                if sum1 > sum2:
                    formed_adn.append(sum1-sum2)
                elif sum1 < sum2:
                    formed_adn.append(sum2-sum1)
            else: continue

        return formed_adn

    def M_nucléosomes():
        #revoie chaque partie de histones pour l'ajouté a un chainon
        le_chainon = []

        for nucléosomes_ind in range(random.choice([2, 4, 6, 8, 10, 12])):
            histones_list_A = nucléosomes.histones()
            histones_list_B = nucléosomes.histones()
            histones_list_C = nucléosomes.histones()
            histones_list_D = nucléosomes.histones()
            histones_list_E = nucléosomes.histones()
            histones_list_F = nucléosomes.histones()
            histones_list_G = nucléosomes.histones()
            histones_list_H = nucléosomes.histones()

            le_chainon.append([histones_list_A, histones_list_B, histones_list_C, histones_list_D, histones_list_E, histones_list_F, histones_list_G, histones_list_H])

        hisone_chainon = []
        for histones in range(len(le_chainon)):
            for histone in range(len(le_chainon[histones])):
                chainon = ""
                for HOPN in range(len(le_chainon[histones][histone])):
                    try:
                        histone_A = le_chainon[histones][histone][HOPN]
                        histone_B = le_chainon[histones][histone+1][HOPN]
                        if histone_A == histone_B:
                            chainon += "A"
                        elif int(histone_A + histone_B) <= 1000:
                            chainon += "B"
                        elif int(histone_A + histone_B) <= 2000:
                            chainon += "C"
                        elif int(histone_A + histone_B) <= 4000:
                            chainon += "D"
                        elif int(histone_A + histone_B) <= 8000:
                            chainon += "E"
                        elif int(histone_A + histone_B) <= 16000:
                            chainon += "F"
                        elif int(histone_A + histone_B) <= 32000:
                            chainon += "G"
                        elif int(histone_A + histone_B) <= 64000:
                            chainon += "H"
                        else:
                            chainon += "Z"
                    except: pass
                if chainon != '':
                    hisone_chainon.append(chainon)

        return hisone_chainon

class adn:
    def interact_A_to_T(adénine, adénine_M, thymine, thymine_M):
        #connect adénine et thymine via deux liaisons hydrogene
        A_T = []
        for i in range(len(adénine_M)):
            A = tool.get_element(adénine_M[i], adénine)
            T =  tool.get_element(thymine_M[i], thymine)
            A_T.append([A, T])
        return A_T

    def interact_G_to_C(guanine, guanine_M, cytosine, cytosine_M):
        #connect guanine et cytosine via trois liaisons hydrogene
        G_C = []
        for i in range(len(guanine_M)):
            G = tool.get_element(guanine_M[i], guanine)
            C = tool.get_element(cytosine_M[i], cytosine)
            G_C.append([G, C])
        return G_C

    def sous_sub_adn():
        #create adn with AT and GC
        adénine = []
        thymine = []
        guanine = []
        cytosine = []

        Ap = adn_part.adénine_A()[0]
        Tp = adn_part.thymine_T()[0]
        Gp = adn_part.guanine_G()[0]
        Cp = adn_part.cytosine_C()[0]

        A_H = tool.indexing_values(Ap, "H")
        A_P = tool.indexing_values(Ap, "P")
        A_O = tool.indexing_values(Ap, "O")
        A_N = tool.indexing_values(Ap, "N")
        #+
        T_H = tool.indexing_values(Tp, "H")
        T_P = tool.indexing_values(Tp, "P")
        T_O = tool.indexing_values(Tp, "O")
        T_N = tool.indexing_values(Tp, "N")

        G_H = tool.indexing_values(Gp, "H")
        G_P = tool.indexing_values(Gp, "P")
        G_O = tool.indexing_values(Gp, "O")
        G_N = tool.indexing_values(Gp, "N")
        #+
        C_H = tool.indexing_values(Cp, "H")
        C_P = tool.indexing_values(Cp, "P")
        C_O = tool.indexing_values(Cp, "O")
        C_N = tool.indexing_values(Cp, "N")

        adénine.append(Ap)
        adénine.append([A_H, A_P, A_O, A_N])
        thymine.append(Tp)
        thymine.append([T_H, T_P, T_O, T_N])
        guanine.append(Gp)
        guanine.append([G_H, G_P, G_O, G_N])
        cytosine.append(Cp)
        cytosine.append([C_H, C_P, C_O, C_N])

        return [adénine, thymine, guanine, cytosine]

    def sub_adn():
        adénine = adn.sous_sub_adn()[0][0]
        thymine = adn.sous_sub_adn()[1][0]
        guanine = adn.sous_sub_adn()[2][0]
        cytosine = adn.sous_sub_adn()[3][0]

        adénine_M = adn.sous_sub_adn()[0][1]
        thymine_M = adn.sous_sub_adn()[1][1]
        guanine_M = adn.sous_sub_adn()[2][1]
        cytosine_M = adn.sous_sub_adn()[3][1]

        A_T = adn.interact_A_to_T(adénine, adénine_M, thymine, thymine_M)
        G_C = adn.interact_G_to_C(guanine, guanine_M, cytosine, cytosine_M)

        A_T_TAB = []
        A_T_TAB.append(A_T)

        G_C_Tab = []
        G_C_Tab.append(G_C)

        return [A_T_TAB, G_C_Tab]

    def get_acide_HPON():
        molecule_H = tool.indexing_values(adn.sub_adn(), "H")
        molecule_P = tool.indexing_values(adn.sub_adn(), "P")
        molecule_O = tool.indexing_values(adn.sub_adn(), "O")
        molecule_N = tool.indexing_values(adn.sub_adn(), "N")

        return [molecule_H, molecule_P, molecule_O, molecule_N]

    def M_adn(number_of_adn: int, molecule_H, molecule_P, molecule_O, molecule_N):
        print(f"[ADN] number of cell for division is {int(number_of_adn)}")
        counting = int(number_of_adn/8)

        complet_adn = []
        for adn_create in range(int(number_of_adn)):
            adn_formé = []
            for org in range(random.choice([8, 10, 12, 14, 16])):
                molecule_adn = []
                for x in range(random.choice([4, 6, 8])):
                    brain_adn = [molecule_H, molecule_P, molecule_O, molecule_N]
                    molecule_adn.append(brain_adn)
                adn_formé.append(molecule_adn)

            full_molecule_branch = []
            for part_molecule in adn_formé:
                int_list = []
                for parted_part_of_molecule in part_molecule:
                    for HPON in parted_part_of_molecule:
                        int_list.append(HPON)
                full_molecule_branch.append(int(sum(int_list)/10))

            if len(full_molecule_branch) >= 15:
                if not full_molecule_branch in complet_adn:
                    complet_adn.append(full_molecule_branch)
                    #print(full_molecule_branch)
                elif counting >= 1000:
                    break
                else:
                    counting += 1

        return complet_adn

class adn_part:
    def adénine_A():
        #generé une modelcule d'adénine (A)
        #une adénine interagit avec une thymine à travers deux liaisons hydrogène
        adénine = [["H", "O"], ["O", "P"], "O", "O", ["O", "H"], "O", [["N", "N"], ["N", "N"]], ["H", "N", "N"]]
        H_lisaison = tool.indexing_values(adénine, "H")
        return [adénine, H_lisaison]

    def cytosine_C():
        #generé une modelcule de cytosine (C)
        cytosine = [["H", "O"], ["O", "P"], "O", "O", ["O", "H"], "O", [[["N", "N"], "O"], [["N", "H"], ["N", "H"]]]]
        H_lisaison = tool.indexing_values(cytosine, "H")
        return [cytosine, H_lisaison]

    def guanine_G():
        #generé une modelcule de guanine (G)
        #une guanine interagit avec une cytosine à travers trois liaisons hydrogène
        guanine = [["H", "O"], ["O", "P"], "O", "O", ["O", "H"], "O", [["N", "N"], ["N", ["N", "H", "N", "H"], "N", "H"], "O"]]
        H_lisaison = tool.indexing_values(guanine, "H")
        return [guanine, H_lisaison]

    def thymine_T():
        #generé une modelcule de thymine (T)
        thymine = [["H", "O"], ["O", "P"], "O", "O", ["O", "H"], "O", [["N", ["O", ["N", "H"]], "O"]]]
        H_lisaison = tool.indexing_values(thymine, "H")
        return [thymine, H_lisaison]

for cell_souche in range(10):
    celulle_count += 1
    celulle_ = cellule.M_cellule()
    print(celulle_)
    with open(f"cellule/cell{celulle_count}.txt", "w") as cell:
        for cell_ in celulle_:
            cell_fis = str(cell_).replace("'", "").replace('[', '').replace(']', '')
            new_cells = cell_fis.split(", ")
            for cells_ in new_cells:
                cell.write(f"{cells_}\n")