from sage.libs.pari.convert_sage import gen_to_sage
import time
import math
import random


# Variable tampon globale pour les polynômes modulaires.
_polmod={}

# Renvoie le polynôme modulaire de degré d.
"""def polmod(d,K):
    A.<X,Y>=PolynomialRing(K,2)
    if not d in _polmod.keys():
        a = pari("polmodular("+str(d)+")")
        #convertir en polynome en 2 variables de sage
        b=gen_to_sage(a, {'x':Y, 'y':Z})
        _polmod[d]=b
    return _polmod[d]"""

def polmod(degres,K):
    """retourne un dictionnaire contenant les polynomes modualaires sur le corps K
    degres= la liste des degres des polynomes que l'on veut"""
    A.<X,Y>=PolynomialRing(K,2)
    P={}
    for d in degres:
        r = pari("polmodular("+str(d)+")")
        s = gen_to_sage(r, {'x':X, 'y':Y})
        P[d] = s
    return P

def get_supersing(k):
    """renvoie une courbe elliptique supersingulière définie sur le corps k"""
    while (True):
        a=k.random_element()
        b=k.random_element()
        if 4*a^3 + 27*b^2 != 0:
            E=EllipticCurve(k, [a, b])
            if E.is_supersingular():
                return E


def voisins(P,sommet,degres=None):
    """P=dictionnaire de polynomes modulaires
    renvoie les voisins de sommet dans le graphe des isogénies de degré donné i.e donne les racines de p(X,sommet)
    pour p dans le dictionnaire p
    la fonction stocke également le degré du lien entre sommet et ses voisins dans le dictionnaire optionnel degres"""
    r=[]
    e=not degres==None
    for d in P:
        X=P[d].parent().gens()[1]
        u=P[d](sommet,X).univariate_polynomial()
        for jm in u.roots():
            if jm!=sommet:
                j=jm[0]
                r.append(j)
            if e:
                degres[(sommet,j)]=d
                degres[(j,sommet)]=d
    return r

def recherche_chemin_depuis_origine(dico,point):
    """donne un chemin allant de l'origine du graphe contenu dans le dictionnaire dico vers point"""
    chemin=[]
    t=point
    while t in dico.keys():
        chemin.append(t)
        t=dico[t]
    chemin.append(t)
    chemin.reverse()
    return chemin

def recherche_cycle(dico,som,voisin):
    """renvoie un cycle dans le graphe contenu dans le dictionnaire dico
    ce cycle démmarre forcémenent à l'origine du parcours et passe par voisin"""
    chemin1=recherche_chemin_depuis_origine(dico,som)
    chemin2=recherche_chemin_depuis_origine(dico,voisin)
    chemin2.reverse()
    fin=len(chemin2)-1
    return chemin1 + chemin2

def liste_degres_depuis_dico(liste,dico):
    """retourne une liste de degrés reliant les j invariants passés dans la liste"""
    degres_liste=[]
    for i in range(0,len(liste)-1):
        degres_liste.append(dico[(liste[i],liste[i+1])])
    return degres_liste


def cycle_parcours_largeur(liste_degres,depart):
    """depart= un j-invariant d'une classe d'iso de courbes supersingulière
    liste_degres= liste des degres qui va constituer le graphe d'isogénies
    La fonction renvoie  un cycle dans le graphe avec pour origine le j-invariant depart"""

    k=depart.parent()
    P=polmod(liste_degres,k)
    dejavu=[]
    a_voir=[depart]
    predecesseur={}
    degres={}
    while (len(dejavu)<k.cardinality() and a_voir):
        sommet=a_voir.pop()
        dejavu.append(sommet)
        successeurs=voisins(P,sommet,degres)
        random.shuffle(successeurs)
        for vois in successeurs:
            if vois in dejavu:
                if vois!=predecesseur[sommet]:
                    res= recherche_cycle(predecesseur,sommet,vois)
                    list_deg=liste_degres_depuis_dico(res,degres)
                    return [res,list_deg]
            else:
                predecesseur[vois]=sommet
                a_voir.append(vois)    
    assert False

def vers_fp(liste_degres,depart):
    """ depart= un j-invariant d'une classe d'iso de courbes supersingulière définie sur Fp^2
    La fonction renvoie un chemin de depart vers un j-invariant défini sur Fp, ce chemin est pris dans 
    le graphe des isogénies de degré liste_degres si celui ci est de taille inférieure à 2*logp et renvoie 0 sinon
    elle renverra également la liste des degrés reliant les j-invariants """

    k=depart.parent()
    Fp=depart.parent().subfield(1)
    logp=math.log2(Fp.cardinality())
    P=polmod(liste_degres,k)
    a_voir=[]
    dejavu=[]
    a_voir.append(depart)
    degres={}
    predecesseur={}
    i=0
    while (a_voir): 
        if(i>=2*logp):
            break
        sommet=a_voir.pop()
        dejavu.append(sommet)
        successeurs=voisins(P,sommet,degres)
        random.shuffle(successeurs)  #ici on melange successeurs pour éviter d'avoir toujours le meme long chemin
        for vois in successeurs:
            if vois in Fp :
                chemin=recherche_chemin_depuis_origine(predecesseur,sommet)
                chemin.append(vois)
                list_deg=liste_degres_depuis_dico(chemin,degres)
                return [chemin,list_deg]
            elif vois not in dejavu:
                a_voir.append(vois)
                predecesseur[vois]=sommet
        i=i+1
    return 0

def vers_fp_court(liste_degres,depart):
    """utilise la fonction vers_fp pour renvoyer un chemin de taille < 2*log(p) allant de depart vers une courbe sur Fp
    elle renverra également la liste des degrés reliant les j-invariants"""

    k=depart.parent()
    i=0
    while(i<k.cardinality()):
        liste=vers_fp(liste_degres,depart)
        if liste!=0:
            return liste
        i=i+1
    assert False


def cycliquedepuisfp(liste_degres,depart): 
    """renvoie un cycle dans le graphe des isogénies de degré liste_degres et partant de depart,
    le procédé de la fonction est de trouver un chemin depart -> j1 -> j2 -> ... -> jn où jn est dans Fp
    on trouve ensuite un chemin   depart -> j'1 -> j'2 -> ... -> j'm où j'm est dans  Fp
    On renvoie alors le chemin depart -> j1 -> j2 -> ... -> jn -> (j(n-1))^p -> ... -> (j1)^p -> depart^p -> (j'(1))^p -> ...
    -> (j'(m-1))^p -> j'm -> j'(m-1) -> ... -> j'1 -> depart"""

    a=vers_fp_court(liste_degres,depart)
    b=vers_fp_court(liste_degres,depart)
    
    liste1=a[0]
    liste2=b[0]
    cycle1=[]
    cycle2=[]
    lon1=len(liste1)
    lon2=len(liste2)
    K=depart.parent()
    p=sqrt(K.cardinality())

    for j in liste1:
        cycle1.append(j)
    for i in range (1,lon1):
        cycle1.append(K(liste1[lon1-i-1]^p))
        """a[1][(cycle1[len(cycle1)-1],cycle1[len(cycle1)-2])]=a[1][(liste1[lon1-i],liste1[lon1-i-1])]
        a[1][(cycle1[len(cycle1)-2],cycle1[len(cycle1)-1],)]=a[1][(liste1[lon1-i],liste1[lon1-i-1])]"""
        a[1].append(a[1][lon1-i-1])
    
    for j in liste2:
        cycle2.append(j)
    for i in range (1,lon2):
        cycle2.append(K(liste2[lon2-i-1]^p))
        """b[1][(cycle2[len(cycle2)-1],cycle2[len(cycle2)-2])]=b[1][(liste2[lon2-i],liste2[lon2-i-1])]
        b[1][(cycle2[len(cycle2)-2],cycle2[len(cycle2)-1])]=b[1][(liste2[lon2-i],liste2[lon2-i-1])]"""
        b[1].append(b[1][lon2-i-1])


    cycle2.reverse()
    return [cycle1 + cycle2[1:], a[1]+b[1]]

def cycle(liste_degres,depart):
    """renvoie un cycle dans le graphe des isogénies de degré liste_degres et partant de depart, la méthode pour trouver ce cycle
    est choisie arbitrairement entre la méthode en passant par fp et celle du parcours du graphe en largeur"""
    choix=random.randint(0,1)
    if choix==0:
        return cycle_parcours_largeur(liste_degres,depart)
    return cycliquedepuisfp(liste_degres,depart)


def calcule_isogenie(E1,j,degre):
    """calcule une isogénie entre la courbe E1 et une courbe E2 dont le j-invariant est donné par j en paramtètre, le degré est en
    paramètre"""
    assert degre <= 3
    divpol=E1.division_polynomial(degre)
    for k in divpol.factor():
        if k[0].degree()==1: #toujours vrai dans le cas supersingulier sur fp2 et degre=2 ou 3?
            phi=E1.isogeny(k[0])
            E2=phi.codomain()
            if E2.j_invariant()==j:
                return phi
    assert False

def isogenie_complete(courbes,degres_entre_courbes,corps): #peut-être remplacer la liste de degrés par le dico des degrés
    """retourne une liste d'isogénie entre des courbes données en argument et sur un corps donné,
    si la liste de courbes donnée est un cycle, la fonction ajoute en fin de liste un isomorphisme entre la courbe de départ et
    celle d'arrivée qui a été donnée par la fonction calcule_isogenie
    ATTENTION, pour calculer une isogénie il a fallu créer une courbe via la méthode EllipticCurve_from_j et donc 
    il faut bien voir que si l'on veut faire des calculs avec une autre courbe de même j-invariant il faudra repasser par un isomorphisme"""
    iso=[]
    E=EllipticCurve_from_j(courbes[0])
    E0=E.change_ring(corps)
    i=0
    for j in courbes[1:]:
        phi=calcule_isogenie(E0,j,degres_entre_courbes[i])
        E0=phi.codomain()
        iso.append(phi)
        i=i+1
    if iso[len(iso)-1].codomain().j_invariant()==iso[0].domain().j_invariant():
        E1=iso[0].domain()
        E2=iso[len(iso)-1].codomain()
        isomo=E2.isomorphism_to(E1)
        iso.append(isomo)
    
    return iso

def evaluer_liste_isogenie(isogenies,point):
    """renvoie l'image d'un point appartenant à une courbe E1 par une liste d'isogénies E1->E2->...->En"""
    for iso in isogenies:
        point=iso(point)
    return point


def elem_ordre_p(E,p,corps):
    """renvoie un élément d'ordre p de la courbe E définie sur le corps donné en argument"""
    E2=E.change_ring(corps)
    neutre=E2(0)
    points=neutre.division_points(p)
    for point in points:
        if point!=neutre:
            return point

    
def duale_de_liste(liste):
    """renvoie l'isogénie duale sous forme de liste d'une isogénie elle-même donnée par une liste en argument"""
    duale=[]
    for i in liste:
        j=i.dual()
        duale.append(j)
    duale.reverse()
    return duale

def degredeliste(liste):
    """renvoie le degré d'une isogénie donnée par une liste d'isogénies en argument"""
    prod=1
    for iso in liste:
        prod=prod*iso.degree()
    return prod

def replace_iso_audessus(isogenies,corps):
    """replace une liste d'isogénies définies entre des courbes sur un certain corps inconnu vers les mêmes courbes définies sur un
    extension"""
    nouvelles=[]
    for iso in isogenies:
        E=iso.domain()
        E1=E.change_ring(corps)
        kerpol=iso.kernel_polynomial()
        pol2=kerpol.change_ring(corps)
        phi=E1.isogeny(pol2)
        nouvelles.append(phi)
    return nouvelles

def produit_des_elements_de_liste(liste):
    prod=1
    for elem in liste:
        prod=prod*elem
    return prod


def produit_scalaire(E,cycle_alpha,cycle_beta,degres_alpha,degres_beta,deg={},ext={}):
    """renvoie le produit scalaire de deux endomorphismes alpha et beta de la courbe E,
    cycle_alpha et cycle_beta désignent la liste des j_invariants à parcourir pour faire un endomorphisme de E 
    de manière déterministe avec la méthode isogenie_complete, il faut donc donner les degrés en argument en plus"""
    # Optional arguments are:
    # - deg, a dictionary of primes l, indexed by the degree of the smallest extension where a point of order l lies
    # - ext, a dictionary of extensions of the base field, indexed by their degree
    K=cycle_alpha[0].parent()
    liste_premiers=[]
    prod=1
    last=1
    resoudre=[]
    torsions=[]
    i=0 
    degalpha=produit_des_elements_de_liste(degres_alpha)
    degbeta=produit_des_elements_de_liste(degres_beta)

    print("NEW",degalpha*degbeta)
    while(prod<=degalpha*degbeta):
        if deg:
            d=min(deg.keys())
            last=deg[d].pop()
            if not deg[d]: del deg[d]
        else:
            last=last.next_prime()
            d=min([x[0].degree() for x in E.division_polynomial(last).factor()])

        if not 2*d in ext:
            ext[2*d]=K.extension(2*d)

        K2=ext[2*d]
        P=elem_ordre_p(E,last,K2)

        alphatemp=isogenie_complete(cycle_alpha,degres_alpha,K2)
        betatemp=isogenie_complete(cycle_beta,degres_beta,K2)
        alpha_du=duale_de_liste(alphatemp)
        beta_du=duale_de_liste(betatemp)
        Qab=evaluer_liste_isogenie(alphatemp,evaluer_liste_isogenie(beta_du,P))
        Qba=evaluer_liste_isogenie(betatemp,evaluer_liste_isogenie(alpha_du,P))
        Q=Qab+Qba
        resoudre.append(discrete_log(Q,P,last,operation='+'))
        prod=prod*last
        liste_premiers.append(last)
        i=i+1 
    return CRT_list(resoudre,liste_premiers)
   


def sous_liste_random(liste):
    """renvoie une sous-liste de taille aléatoire de la liste en argument"""
    copie=deepcopy(liste)
    res=[]
    taille=randint(1,len(liste))
    i=0
    while i<taille and copie:
        random.shuffle(copie)
        elem=copie.pop()
        res.append(elem)
        i=i+1
    return res
        


def matrice_gram(E):
    """renvoie la matrice de Gram d'une courbe E ainsi que les endomorphismes qui servent à la construire"""
    ext={}
    deg={}
    for l in primes(62):
        d=min([x[0].degree() for x in E.division_polynomial(l).factor()])
        if not d in deg.keys(): deg[d]=[]
        deg[d].append(l)

    MS=MatrixSpace(ZZ,9)
    jinv=E.j_invariant()
    K=jinv.parent()

    base=[]
    Gram=MS(0)
    lis=[2,3]
    while Gram.rank()<4 and len(base)<10:
        
        souslis=sous_liste_random(lis)
        a=cycle(souslis,jinv)
        cpt=0
        while produit_des_elements_de_liste(a[1])>=31381059609 and cpt<10: #Borne arbitraire pour eviter une execution trop longue sur les prod scal
            a=cycle(souslis,jinv)
        base.append(a)
        for i in range(len(base)):
            cur=len(base)-1
            newdeg=deepcopy(deg)
            Gram[i,cur]=produit_scalaire(E,base[i][0],base[cur][0],base[i][1],base[cur][1],deg=newdeg,ext=ext)
            Gram[cur,i]=Gram[i,cur]
    assert Gram.rank()==4

    return [Gram,base]

def base_depuis_matrice_gram(matrice,elements):
    """renvoie la forme normale de hermite de la matrice de Gram et une base de l'anneau d'endomorphismes
    matrice=matrice de Gram
    elements= endomorphismes qui servent à construire la matrice de Gram"""

    #U=matrice.LLL_gram()
    #res=U.transpose() * matrice * U
    Hermite=matrice.hermite_form()
    base=[]
    for i in range(len(elements)):
        E=0
        for j in range(len(elements)):
            E= E+Hermite[j][i]*elements[i] #on considère que les colonnes de la matrice engendre le sous-réseau de l'anneau d'endomorphismes
        base.append(E)
    return [Hermite,base]


def main():
    p=863
    K.<c>=GF(p^2,name='c',modulus=x^2 -x +5)

    P=polmod([2,3],K)
    
    #E1=EllipticCurve(K,[40*c+535,720*c+768])
    E=EllipticCurve_from_j(K(403*c+428))
    print(E.j_invariant(),E)
    assert E.is_supersingular()
    #w=cycle_parcours_largeur([2,3],E.j_invariant())
    #print(w[0])
    #print(w[1])
    #isoge=isogenie_complete(w[0],w[1],K)
    #print(isoge)
    #a=produit_scalaire(E,w[0],w[0],w[1],w[1])
    #print(a,2*produit_des_elements_de_liste(w[1]))
    
    MS=MatrixSpace(ZZ,9)
    
    
    M=matrice_gram(E)
    print(M)
    print(base_depuis_matrice_gram(M[0],M[1]))




if __name__ == '__main__':
    main()