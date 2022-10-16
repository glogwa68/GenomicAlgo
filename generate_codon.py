def trad_ARN(transcrit):
    traduit = ""

    # **** Boucle mystère : à quoi sert-elle ? ****
    initiation = False
    numbis = 0
    while numbis < len(transcrit)-2:
        if transcrit[numbis:numbis+3] == "AUG":
            initiation = True
            break
        numbis+= 1
    # **** Fin de la boucle mystère ****

    # Code génétique
    ala = ["GCU", "GCC", "GCA", "GCG"]
    arg = ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"]
    asn = ["AAU", "AAC"]
    asp = ["GAU", "GAC"]
    cys = ["UGU", "UGC"]
    gln = ["CAA", "CAG"]
    glu = ["GAA", "GAG"]
    gly = ["GGU", "GGC", "GGA", "GGG"]
    his = ["CAU", "CAC"]
    ile = ["AUU", "AUC", "AUA"]
    leu = ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"]
    lys = ["AAA", "AAG"]
    met = ["AUG"]
    phe = ["UUU", "UUC"]
    pro = ["CCU", "CCC", "CCA", "CCG"]
    ser = ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"]
    thr = ["ACU", "ACC", "ACA", "ACG"]
    trp = ["UGG"]
    tyr = ["UAU", "UAC"]
    val = ["GUU", "GUC", "GUA", "GUG"]
    stop = ["UAG", "UGA", "UAA"]

    # **** Boucle pour traduire la séquence ****
    if initiation:
        while numbis < len(transcrit)-2:
            if transcrit[numbis:numbis+3] in ala:
                traduit += "Ala"
            elif transcrit[numbis:numbis+3] in arg:
                traduit += "Arg"
            elif transcrit[numbis:numbis+3] in asn:
                traduit += "Asn"
            elif transcrit[numbis:numbis+3] in cys:
                traduit += "Cys"
            elif transcrit[numbis:numbis+3] in gln:
                traduit += "Gln"
            elif transcrit[numbis:numbis+3] in glu:
                traduit += "Glu"
            elif transcrit[numbis:numbis+3] in gly:
                traduit += "Gly"
            elif transcrit[numbis:numbis+3] in his:
                traduit += "His"
            elif transcrit[numbis:numbis+3] in ile:
                traduit += "Ile"
            elif transcrit[numbis:numbis+3] in leu:
                traduit += "Leu"
            elif transcrit[numbis:numbis+3] in lys:
                traduit += "Lys"
            elif transcrit[numbis:numbis+3] in met:
                traduit += "Met"
            elif transcrit[numbis:numbis+3] in phe:
                traduit += "Phe"
            elif transcrit[numbis:numbis+3] in pro:
                traduit += "Pro"
            elif transcrit[numbis:numbis+3] in ser:
                traduit += "Ser"
            elif transcrit[numbis:numbis+3] in thr:
                traduit += "Thr"
            elif transcrit[numbis:numbis+3] in trp:
                traduit += "Trp"
            elif transcrit[numbis:numbis+3] in tyr:
                traduit += "Tyr"
            elif transcrit[numbis:numbis+3] in val:
                traduit += "Val"
            elif transcrit[numbis:numbis+3] in stop:
                break
            numbis += 3

        print("Séquence traduite :", traduit)
        return traduit
    else:
        return False