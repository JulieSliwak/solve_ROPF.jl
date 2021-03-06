# solve_ROPF.jl

## Problème OPF avec activation possible des MCS et contrainte sur la production active (contraintes GENmoves)
2 fichiers nécessaires :
- une instance MATPOWER (.m)
- un plan de production active (fichier .csv avec 5 colonnes (nom du générateur ; P ; ; Pmin ; Pmax) par exemple (Gen_Sgen_183_Re ; 8.98 ; ; 7.98 ; 26.6)

Un fichier .dat est construit à partir de l'instance MATPOWER .m dans la fonction `construct_dat_file_ROPF(instance_name, matpower_instance_path, output_path)` dont les arguments sont :
* `instance_name` : le nom de l'instance MATPOWER, par exemple "case14"
* `matpower_instance_path` : le chemin de l'instance MATPOWER par exemple "D:\\repo\\data\\matpower\\case14.m"
* `output_path` : le chemin du répertoire dans lequel stocker les instances .dat par exemple "D:\\repo\\data\\data_ROPF\\RTE_ROPF"

La structure `ROPF_infos` a été définie pour résumer toutes les informations à propos d'un problème ROPF :
* `instance_name` : le nom de l'instance MATPOWER, par exemple "case14"
* `matpower_instance_path` : le chemin de l'instance MATPOWER par exemple "D:\\repo\\data\\matpower\\case14.m"
* `output_path` : le chemin du répertoire dans lequel sont stockées les instances .dat par exemple "D:\\repo\\data\\data_ROPF\\RTE_ROPF"
* `decomposition` : le nom de la décomposition en cliques maximales à utiliser, par exemple "cholesky"
* `output_decomposition_path` : le chemin du répertoire dans lequel sont stockées les décompositions par exemple "D:\\repo\\data\\data_sdp"
* `generation` : le nom du plan de production active à utiliser, par exemple "generationrandom1"
* `generation_files_path` : le chemin du répertoire dans lequel sont stockés les plans de production par exemple "D:\\repo\\data\\data_ROPF\\RTE_ROPFu"

## Relaxation SDP et décomposition en cliques
Fonctions dans le fichier `solve_SDP.jl`

Si aucune décomposition en cliques maximales n'est connue, une décomposition en cliques maximales est calculée avec une factorisation de Cholesky précédée d'un ordre AMD avec la fonction `generate_clique_decomposition(instance_name, matpower_instance_path, output_decomposition_path)` dont les arguments sont :
* `instance_name` : le nom de l'instance MATPOWER, par exemple "case14"
* `matpower_instance_path` : le chemin de l'instance MATPOWER par exemple "D:\\repo\\data\\matpower\\case14.m"
* `output_decomposition_path` : le chemin du répertoire dans lequel stocker les décompositions par exemple "D:\\repo\\data\\data_sdp"
La décomposition calculée s'appelle alors "cholesky".

La fonction `solve_SDP(ROPF, flag)` résout la relaxation SDP du problème `ROPF` qui a une structure `ROPF_infos` soit en version montée du plan de production (si `flag="plus"`) soit en version baisse du plan de production (si `flag="minus"`). Le log est enregistré dans le répertoire "Mosek_runs".


## Encadrement de la valeur optimale
La fonction `UB_plus, LB_plus, UB_minus, LB_minus = solve1(ROPF)` retourne un encadrement de la solution du problème `ROPF` (borne supérieure (UB) et borne inférieure (LB) pour chaque version du problème (baisse ou montée du plan de production)).

Le problème `ROPF` est traité en deux temps : on étudie la version montée du plan de production (`flag="plus"`) et la version baisse du plan de production (`flag="minus"`).

On commence par calculer une borne inférieure pour chacune des deux versions en résolvant la relaxation SDP.

Si l'une des relaxations SDP est non réalisable, cela implique que le problème est non réalisable donc la borne supérieure retournée est `UB=Inf`. Sinon, une solution réalisable est calculée en trois étapes avec Knitro via AMPL :
* résolution de la relaxation continue du problème pour déterminer un bon point de départ pour la résolution du problème avec les variables binaires
* résolution du problème avec les variables binaires en utilisant l'option MPEC de Knitro (contraintes de complémentarité)
* fixation des variables binaires en les arrondissant puis résolution du problème continu



## Algorithme B&B
La fonction `UB_plus, LB_plus, UB_minus, LB_minus = solve2(ROPF, max_time)` retourne un encadrement potentiellement plus précis de la solution du problème `ROPF` car une procédure de Branch-and-Bound (B&B) est utilisée pour calculer une solution réalisable. Cette procédure
est basée sur une relaxation SDP et ne concerne que les variables binaires (ce n'est pas un B&B spatial donc pas une méthode exacte). Pour éviter que le B&B ne soit pas trop long, les variables binaires au-dessus et en-dessous d'un certain seuil sont fixées (< 10^(-4) et > 0.9 par défaut) et une limite de temps doit être précisée (`max_time`).

## Exemple
Fichier `solve_Matpower_instance.jl`
* Si l'instance MATPOWER n'a jamais été utilisée avant, utiliser la fonction `construct_dat_file_ROPF(instance_name, matpower_instance_path, output_instance_path)` pour construire l'instance OPF en format .dat (construction du problème en variables complexes à partir du fichier MATPOWER, conversion en problème réel avec la représentation cartésienne, export en fichier texte). Cette étape peut être longue, en particulier la fonction de conversion est lente dès que le problème contient plus d'une centaine de noeuds.
* Si aucune décomposition en cliques n'est connue, utiliser la fonction `generate_clique_decomposition(instance_name, matpower_instance_path, output_decomposition_path)` pour calculer une décomposition à partir d'une factorisation de Cholesky précédé d'un ordre AMD.
* Le plan de production active doit avoir été généré précédemment (avec la bonne syntaxe pour le nom des générateurs : e.g. "Gen_Sgen_2_Re" pour le générateur au noeud 2).
* Une fois que le fichier .dat existe et qu'il existe aussi une décomposition en cliques maximales, utiliser la fonction `UB_plus, LB_plus, UB_minus, LB_minus = solve1(ROPF)` pour obtenir un premier encadrement de la valeur optimale et/ou la fonction `UB_plus, LB_plus, UB_minus, LB_minus = solve2(ROPF, max_time)` pour obtenir un encadrement potentiellement plus précis.




## Packages nécessaires
Julia Version 1.0.3
* Mosek.jl avec Mosek installé
* MathProgComplex.jl (code personnel, branche dev 1.0.3)
* ComplexOPF.jl (code personnel)
* LightGraphs.jl
* MetaGraphs.jl
* LinearAlgebra.jl
* SparseArrays.jl
* SuiteSparse.jl
* DelimitedFiles.jl
* KNITRO.jl
* Xpress.jl
* JuMP.jl
* DataStructures.jl
* Printf.jl
* ArgPars.jl

* AMPL, KNITRO et Xpress


## Instructions
Clôner les trois dépôts (idéalement au même endroit dans un répertoire que l'on notera `repo`) :
* MathProgComplex.jl (Attention : aller sur la branche dev 1.0.3)
* ComplexOPF.jl
* solve_ROPF_GENmoves.jl

Modifier la première ligne du fichier `repo\ComplexOPF.jl\src\PowSysMod_body.jl`pour mettre à jour le chemin où se trouve le package MathProgComplex.jl. Si le package se trouve dans le même répertoire que les deux autres, la commande `push!(LOAD_PATH, dirname(pwd()))` peut être utilisée. Sinon spécifier le chemin : `push!(LOAD_PATH, path_to_MathProgComplex.jl_package)`.

Aller dans le répertoire solve_ROPF_GENmoves.jl

Modifier la première ligne du fichier `generic_functions.jl` pour mettre à jour le chemin du fichier `\ComplexOPF.jl\src\PowSysMod_body.jl`. Par exemple : `include(joinpath("..", "ComplexOPF.jl","src", "PowSysMod_body.jl"))` si ComplexOPF.jl et solve_ROPF_GENmoves.jl se trouvent dans le même répertoire

Créer un répertoire avec les instances Matpower à traiter, par exemple `data_matpower`.

Créer un répertoire avec les fichiers csv pour la production, par exemple `data_ROPF`.

Modifier le fichier `solve_Matpower_instance.jl` pour donner le nom de l'instance Matpower à traiter et les chemins des différents fichiers à utiliser.
