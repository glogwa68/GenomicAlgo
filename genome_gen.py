import random
import numpy as np
import base64
import string
import threading
import tools
import tools as tool
import queue

global recurent
recurent = 0
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

"""class cellule:
    def membrane(self):
        #gen taille fixe de membrane pour cellule

    def noyau(self):
        #chromosome

    def M_cellule(self):
        #noyau
        #membrane

class chromosome:
    def generate_double_chormatides(self):
        #chromatides x2 (assemblé)

    def M_chromosome(self):
        #generate double chromatides add to tab 2 +100 a chaque telomere


class chromatides:
    def telomere(self):
        #double brins ADN

    def centromere(self):
        #connect fin de brin telomere to exact positon de debut telomere inversé superieur x4

    def M_chromatides(self):
        #generate 2 tab with 4 conection after 1er=> aux mileur a droite, 2eme=> milieux gauche connect 4 branche of 1/2 per tab
"""
class nucléosomes:
    def calcul_adn(division_int: int):
        adn_fullpart = adn.M_adn(int(division_int))  # range_of_adn
        print(len(adn_fullpart))
        part_adn = []
        for molecule in adn_fullpart:
            divided_molecule = tool.diviser_liste(molecule, random.randint(4096, 8198))
            parted = 0
            for part in divided_molecule:
                for i in range(len(part)):
                    parted += int(part[i])
            part_adn.append(parted)

        sum_adn = []
        divided_adn = tools.diviser_liste(part_adn, random.randint(256, 1024))
        for adn_part in divided_adn:
            sum_adn.append(sum(adn_part))

        print(f"(divided by {division_int}) :: resulted list with {len(sum_adn)} rows.")
        return sum_adn

    def histones():
        #enroule une elise d'adn en bobine de 1 lise d'adn de longeur d'un chaine d'adn et les fusione avec chaque brin.
        #prend en entré l'adn pour le repliqué pour donné des nucléosomes
        range_of_adn = random.randint(10000000, 18000000)

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
        pass

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

    def M_adn(number_of_adn: int):
        new_range = number_of_adn // random.randint(1000, 1800)

        print(f"[ADN] number of molecules for is {int(new_range)}")
        counting = int(new_range/8)

        complet_adn = []
        for adn_create in range(int(new_range)):
            adn_formé = []
            for org in range(random.choice([8, 10, 12, 14, 16])):
                molecule_adn = []
                for x in range(random.choice([4, 6, 8])):
                    molecule_H = tool.indexing_values(adn.sub_adn(), "H")
                    molecule_P = tool.indexing_values(adn.sub_adn(), "P")
                    molecule_O = tool.indexing_values(adn.sub_adn(), "O")
                    molecule_N = tool.indexing_values(adn.sub_adn(), "N")
                    brain_adn = [molecule_H, molecule_P, molecule_O, molecule_N]
                    molecule_adn.append(brain_adn)
                adn_formé.append(molecule_adn)

            full_molecule_branch = []
            for part_molecule in adn_formé:
                int_list = []
                for parted_part_of_molecule in part_molecule:
                    for HPON in parted_part_of_molecule:
                        for parted_HPON in HPON:
                            int_list.append(sum(parted_HPON))
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

print(nucléosomes.histones())