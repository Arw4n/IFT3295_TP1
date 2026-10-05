"""IFT3295 - TP1 - Assemblage de fragments (Ex. 2.2 et 2.3).

Ce fichier est A COMPLETER. Regles:
* Completez uniquement le corps des fonctions ci-dessous.
* Ne modifiez ni les noms, ni les signatures, ni les types.
* Vous pouvez utiliser NetworkX pour representer et manipuler le graphe.
* Typez tous vos parametres et retours (types simples: int, str, list, ...).
* Le programme principal est ``main.py`` (fourni): il appelle les fonctions
  de ce module, il ne faut donc rien executer ici directement.
"""

import networkx as nx

from overlap import chevauchement_maximal

SEUIL_DEFAUT: int = 80

# Les noeuds du graphe sont les entiers 0..n-1 (indices des reads); chaque
# arete porte l'attribut entier ``poids`` (= score du chevauchement).


def graphe_chevauchements(
    scores: list[list[int]], seuil: int = SEUIL_DEFAUT
) -> nx.DiGraph:
    """Construit le graphe oriente des chevauchements pertinents.

    Args:
        scores (list[list[int]]): Matrice carree des scores de chevauchement.
        seuil (int): Score minimal requis pour conserver une arete.

    Returns:
        nx.DiGraph: Graphe dont les noeuds sont les indices ``0..n-1``.
        Une arete ``(i, j)`` existe si ``scores[i][j] >= seuil``; son attribut
        ``poids`` vaut le score correspondant.
    """
    
    G = nx.DiGraph()
    G.add_nodes_from([i for i in range(len(scores))]) # création des noeuds du graphe
    for i in range(len(scores)): # Parcours en O(n²)
        for j in range(len(scores[i])):
            if scores[i][j]>=seuil:
                G.add_edge(i, j, poids=scores[i][j])

    return G
    raise NotImplementedError  #TODO


def reduction_transitive(graphe: nx.DiGraph) -> nx.DiGraph:
    """Calcule la reduction transitive du graphe de chevauchement.

    Supprime une arete ``(u, v)`` lorsqu'un autre chemin simple de ``u`` a
    ``v`` subsiste. Pour un cycle de deux noeuds, conserve uniquement l'arete
    de plus grand poids.

    Args:
        graphe (nx.DiGraph): Graphe oriente dont les aretes portent un attribut
            entier ``poids``.

    Returns:
        nx.DiGraph: Copie reduite du graphe d'entree. Le graphe d'entree n'est
        pas modifie.
    """
    GT = graphe.copy()

    Deux_cycles = [(u, v) for u, v in GT.edges if GT.has_edge(v, u) and u < v] # Vérification de 2-cycles
    for i in Deux_cycles:
        if GT.edges[i[0], i[1]]["poids"] >= GT.edges[i[1], i[0]]["poids"]: # Garder le plus lourd
            GT.remove_edge(i[1], i[0])
        else:
            GT.remove_edge(i[0], i[1])
    
    for origine in list(GT.nodes()): # DFS pour voir les chemins redondants
        for dest in list(GT.successors(origine)):
            visite={origine}
            voisins = [i for i in GT.successors(origine) if i!=dest]

            while len(voisins)!=0:
                act = voisins.pop()
                if act not in visite:
                    visite.add(act)
                    if act == dest:
                        GT.remove_edge(origine, dest)
                        break
                    voisins.extend(GT.successors(act))
    return GT
    raise NotImplementedError  # TODOGT.remove_edge(i[1], i[0])


def ordre_assemblage(graphe: nx.DiGraph) -> list[int]:
    """Trouve l'ordre d'assemblage des reads dans le graphe reduit.

    Args:
        graphe (nx.DiGraph): Graphe reduit produit par
            ``reduction_transitive``.

    Returns:
        list[int]: Indices des reads dans un chemin du graphe qui visite
        chaque noeud exactement une fois. L'enonce garantit l'existence de
        ce chemin apres reduction.
    """
    ordre = [] 
    for i in graphe.nodes: # On teste jusqu'à trouver le bon noeud de départ
        aretes = list(nx.dfs_edges(graphe,source=i)) # Affiche l'arbre du graphe selon un DFS

        sources = [tpl[0] for tpl in aretes]
        apparitions = [tpl[0] for tpl in aretes if sources.count(tpl[0])>1] # Vérifie si un sommet est visité plus d'une fois


        if len(aretes)!=len(graphe.nodes())-1: # Si tous les sommets ne sont pas atteints
            continue
        elif len(apparitions)>0:
            continue
        else:
            ordre = [tpl[0] for tpl in aretes]
            ordre.append(aretes[-1][1]) # Ajoute le dernier sommet
            break

    return ordre
    raise NotImplementedError  # TODO

def sequence_finale(reads: list[str], ordre: list[int]) -> tuple[str, list[int]]:
    """Assemble les reads en une sequence de fragment genomique.

    Args:
        reads (list[str]): Sequences des reads a assembler.
        ordre (list[int]): Indices des reads dans l'ordre d'assemblage.

    Returns:
        tuple[str, list[int]]: Tuple contenant la sequence assemblee et les
        longueurs des chevauchements entre les paires de reads consecutifs,
        dans l'ordre.
    """
    read_tot = reads[ordre[0]]
    l_chevauchements = []

    for i in range(len(ordre)-1):
        read1 =reads[ordre[i]]
        read2 = reads[ordre[i+1]]
        chevauchement = chevauchement_maximal(read1, read2)[3]
        l_chevauchements.append(chevauchement)
        read_tot +=read2[chevauchement:]
    
    return (read_tot, l_chevauchements)
    raise NotImplementedError  # TODO
