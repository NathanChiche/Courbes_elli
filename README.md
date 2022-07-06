# Courbes_elli
Anneau d'endomorphismes de courbes elliptiques supersingulières

Le but de ce code est de calculer l'anneau d'endomorphismes de courbes supersingulières.
Il procure quelques fonctions dans ce but.
L'approche est de calculer des endomorphismes au hasard d'une courbe donnée via différentes manières (calcul d'un cycle dans le graphe des isogénies de certains degrés via un parcours en largeur ou en passant allant jusqu'à Fp puis en remontant avec le Frobenius)
Il procure ensuite des méthodes de calculs sur des endomorphismes donnés sous cette forme spéciale (cycle de j-invariants)
Il peut alors calculer la matrice de Gram en effectuant les produits scalaires entre ces endo.

Il reste dans ce code à exprimer l'anneau d'endomorphismes en tant que réseau de dim 4 et non un sous-réseau ce celui-ci.

A noter que les exécutions sont longues, il vaut mieux procéder par étapes pour faire des tests.
Pour des raisons de facilité, les calculs sont toujours faits à partir d'une courbe E0 créée via la méthode Elliptic_curve_from_j.
Pour obtenir l'anneau d'endomorphismes d'une courbe donnée E1 ou toute autre information, il faut donc penser à utiliser l'isomorphisme reliant ces deux courbes de même j invariant.

