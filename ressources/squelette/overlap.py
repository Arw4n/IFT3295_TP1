"""IFT3295 - TP1 - Algorithme de chevauchement de sequences (Q1.5, Ex. 2.1).

Ce fichier est A COMPLETER. Regles:
* Completez uniquement le corps des fonctions ci-dessous.
* Ne modifiez ni les noms, ni les signatures, ni les types, ni les constantes.
* Typez tous vos parametres et retours (types simples: int, str, list, ...).
* Le programme principal est ``main.py`` (fourni): il appelle les fonctions
  de ce module, il ne faut donc rien executer ici directement.
"""

MATCH: int = 4
MISMATCH: int = -4
INDEL: int = -8


def chevauchement_maximal(x: str, y: str) -> tuple[int, str, str, int]:
    """Calcule le meilleur chevauchement ordonne entre deux sequences.

    Args:
        x (str): Premiere sequence a aligner.
        y (str): Deuxieme sequence a aligner.

    Returns:
        tuple[int, str, str, int]: Score maximal, deux lignes alignees et
        longueur du chevauchement.

    Examples:
        >>> chevauchement_maximal("ACCA", "CACGC")
        (8, 'CA', 'CA', 2)
        >>> chevauchement_maximal("CACGC", "ACCA")
        (4, 'ACGC', 'AC-C', 4)
    """
    m = len(x)
    n = len(y)

    V = [[0] * (n + 1) for _ in range(m + 1)]   # premiere ligne et premiere colonne à 0
    P = [[""] * (n + 1) for _ in range(m + 1)]  # directions

    for i in range(1,m+1):
        for j in range(1,n+1):
            if x[i-1] == y[j-1]:
                match = 4 #match
            else:
                match = -4 #mismatch
            diag = V[i-1][j-1] + match
            haut = V[i-1][j] - 8
            gauche = V[i][j-1] - 8
            score_max, direction = max((diag, "diag"), (haut, "haut"),(gauche, "gauche"))
            V[i][j] = score_max
            P[i][j] = direction


    # trouver valeur max dans la derniere ligne et colonne
    
    max_total = -float("inf")
    i_max, j_max = m, n

    for j in range(1, n + 1): 
        if V[m][j] > max_total:
            max_total = V[m][j]
            i_max, j_max = m, j

    # faire le backtracking
    align_x = ""
    align_y = ""


    while i_max > 0 and j_max > 0:
        direction = P[i_max][j_max] 
    
        if direction == "diag":
            align_x = x[i_max - 1] + align_x
            align_y = y[j_max - 1] + align_y
            i_max -= 1
            j_max -= 1
        elif direction == "haut":
            align_x = x[i_max - 1] + align_x
            align_y = "-" + align_y
            i_max -= 1
        elif direction == "gauche":
            align_x = "-" + align_x
            align_y = y[j_max - 1] + align_y
            j_max -= 1
    longueur_chevauchement = len(align_x)
    return max_total, align_x, align_y, longueur_chevauchement


def matrice_chevauchements(reads: list[str]) -> list[list[int]]:
    """Construit la matrice des scores de chevauchement de toutes les paires.

    Args:
        reads (list[str]): Sequences a comparer, par exemple les 20 reads de
            ``reads.fq``.

    Returns:
        list[list[int]]: Matrice carree dont l'element ``[i][j]`` est le score
        maximal de la paire ordonnee ``(reads[i], reads[j])``. La diagonale
        contient des zeros.
    """
    N = len(reads)

    M = [[0] * N for _ in range(N)]

    for i in range(N):
        for j in range(N):
            if i!=j:
                score_max, _, _, _ = chevauchement_maximal(reads[i], reads[j])
                M[i][j] = score_max
            else:
                M[i][j] = 0     # si les reads sont pareils

    return M