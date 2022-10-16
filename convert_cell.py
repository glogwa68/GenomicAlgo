import generate_codon

all_cells = []
for cells in range(1, 11):
    with open(f"cellule/cell{cells}.txt", "r") as cell:
        cell_tab = cell.readlines()

    converted_cell = []
    for ARN in range(len(cell_tab)):
        if ARN != len(cell_tab)-1:
            #ABCDEFGHZ
            new_ARN = str(cell_tab[ARN]).replace("\n", "").replace("B", "U").replace("H", "C").replace("F", "G").replace("D", "AA").replace("E", "AU").replace("F", "UA").replace("G", "UU").replace("H", "CG").replace("Z", "AA")
            converted_cell.append(new_ARN)

    all_cells.append(converted_cell)

traduite_cell = []
for cell in range(len(all_cells)):
    for cell_ in range(len(all_cells[cell])):
        celling = str(all_cells[cell][cell_]).replace("UUUU", "UGAA")[:33]
        trans_cell = generate_codon.trad_ARN(celling)
        if trans_cell != False:
            traduite_cell.append(trans_cell)

with open("trasnfuge_cells/trans_10X_cells.txt", "w") as transfuge:
    for celling_phase in traduite_cell:
        transfuge.write(f"{celling_phase}\n")



