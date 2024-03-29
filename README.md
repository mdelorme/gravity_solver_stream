# Gravitational solver stream

Dépôt github correspondant aux streams de développements de ma chaîne twitch (https://www.twitch.tv/astro_md).

**Note** A partir de la Session #3, on utilise la bibliothèque de portage de performances [Kokkos](https://github.com/kokkos/kokkos). Consultez le lien pour plus d'informations à ce sujet. Notez en revanche qu'il n'est pas nécessaire d'installer Kokkos sur votre système. La bibliothèque est livrée avec le code du stream sous la forme d'un sous-module git.

## Installation et compilation

### Dépendances
Ces instructions ont été testées sur un Ubuntu 20.10 mais devraient fonctionner pour tout système disposant de :
  * Un compilateur C++ compatible C++11. 
  * `git` pour récupérer le repo local
  * `cmake` en version > 3.16

### Obtenir les sources
Pour cloner le repo localement:

```bash
git clone --recursive https://github.com/mdelorme/gravity_solver_stream
```

Si vous n'avez pas utilisé l'option `--recursive` lors du clonage, il faut initialiser et mettre à jour les sous-modules git pour récupérer Kokkos:

```bash
git submodule init
git submodule update
```

### Compilation
La compilation utilise cmake, pour respecter les bonnes pratiques, on commence par créer un dossier pour la compilation

```bash
mkdir build
cd build
```

On initialise le système de build avec la commande `cmake` : 

```bash
cmake ..
```

Afin d'utiliser Kokkos à partir de la session 3, il faut aussi indiquer sur quel backend parallèle on souhaite compiler le code. Lancer `cmake` sans option indiquer à Kokkos de compiler sur le backend `Serial` (sans parallélisme). Pour que le code soit compilé sur le backend OpenMP (parallélisme multi-thread sur CPU), on rajoute l'option `-DKokkos_ENABLE_OPENMP=ON`. Pour que le code soit compilé sur le backend Cuda on utilise l'option `-DKokkos_ENABLE_CUDA=ON`. Si vous connaissez l'architecture sur laquelle vous souhaitez compiler, vous pouvez aussi l'indiquer dans l'option `-DKokkos_ARCH_{NomDeLArchitecture}=ON`. Une fois l'étape de configuration avec `cmake` effectuée, vous pouvez compiler le code en tapant `make -j`.

**Exemple 1: Compilation serial**

```bash
cmake ..
make -j
```

**Exemple 2: Compilation OpenMP sans info d'architecture**

```bash
cmake -DKokkos_ENABLE_OPENMP=ON ..
make -j
```

**Exemple 3: Compilation Cuda sans info d'architecture**

```bash
cmake -DKokkos_ENABLE_CUDA=ON ..
make -j
```

**Exemple 4: Compilation Cuda pour une architecture Turing**

```bash
cmake -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_TURING75=ON ..
make -j
```

Si vous souhaitez voir la liste des options et que vous disposez de `cmake-ncurses` vous pouvez simplement taper `ccmake ..` et éditer vous même la liste des options disponibles.

## Descriptif des sessions
### Session 1 (11/12/2020)
Lien vidéo : https://www.youtube.com/watch?v=dXc-2P8fHfk

Contenu : 
 * Explication du problème à n-corps. 
 * Présentation de la méthode de sommation directe.
 * Implémentation d'un système à 2 corps, Terre-Soleil
 <img src="https://github.com/mdelorme/gravity_solver_stream/blob/main/figs/session1.png" width="300" height="300" />
 
### Session 2 (15/01/2021) 
Lien vidéo : https://www.youtube.com/watch?v=__KIreszA6I
 
Contenu :
 * Nettoyage du code pour une approche objet
 * Extension à un nombre arbitraire de particules
 * Mise en place de particules traceuses 
 * Exemple jouet sur la collision de deux galaxies
 
 <img src="https://github.com/mdelorme/gravity_solver_stream/blob/main/figs/session2.gif" width="300" height="300" />
 
### Session 3 (05/03/2021)
Lien vidéo : https://youtu.be/UNgI8dhXvs0

Contenu :
 * Présentation des concepts généraux du parallélisme
 * Présentation de la librairie Kokkos
 * Réécriture du code en utilisant la librairie Kokkos
 * Tests avec un système de 2000 particules, comparaison OpenMP/Sérial
 * Parallélisme sur GPU

Comparaison des performances

|    Type de run     | Temps |
|   :-----------:    |:-----:|
| Serial             | 81s   |
| OpenMP, 4 threads  | 37s   |
| OpenMP, 8 threads  | 20s   |
| OpenMP, 16 threads | 20s   |
| GPU                | 15s   |

Tests réalisés sur :
  * CPU: Intel(R) Xeon(R) E-2286M  CPU @ 2.40GHz.
  * GPU: Nvidia Quadro T2000

Le code tourne désormais en parallèle sur différents backends mais n'est pas encore optimisé !

### Session 4 (02/04/2021)

Contenu :
  * Optimisation du code pour GPU
  * Introduction de parallélisme hiérarchique
  * Run plus avec un plus grand nombre de particules massives

Résultats finaux : Échec complet d'optimisation. Plus d'investigations sont nécessaires. On fera peut être ça plus une autre fois. Je n'uploade pas la vidéo, mais le contenu des tests est dispo dans le dossier `session4`.


### Session 5 (09/07/2021)

Lien vidéo : https://www.youtube.com/watch?v=jKL9OkwUxKQ

Contenu :
  * Introduction à l'algorithme de Barnes & Hut
  * Codage de l'algorithme, sur CPU
  * Runs divers et variés, et comparaison avec l'algorithme de sommation directe.


#### Comparaison
Résultats finaux des temps de calcul pour 100 itérations, avec OpenMP sur 16 threads :

|Nombre de particules|Sommation directe|Barnes & Hut|
|:------------------:|:---------------:|:----------:|
|100                 | 0.037s          | 0.047s     |
|1000                | 0.301s          | 0.300s     |
|10000               | 14.734s         | 1.818s     |
|100000              | 1324.618s       | 19.699s    |
|1000000             | N/A             | 197.954s   |


#### Représentation graphique :

<img src="https://github.com/mdelorme/gravity_solver_stream/blob/main/figs/session4.png" width="512" height="512" />
En bleu la méthode de sommation directe, en rouge la méthode de Barnes & Hut. En pointillés, les projections des complexités théoriques des deux algorithmes. On voit que très rapidement, les deux courbes s'alignent sur la tendance théorique. En projetant, on obtient donc un temps de calcul de ~10^5 secondes pour faire 100 itérations à 1 million de points en sommation directe, ce qui correspond environ à 27h de calcul.


#### Film, pour 10K particules :

<img src="https://github.com/mdelorme/gravity_solver_stream/blob/main/figs/session4.gif" width="400" height="400" />

## Session 6 (20/10/2021) :
 
Contenu :
  * Reprise douce et progressive
  * Gestion des paramètres de manière plus propre
  * Gestion des unités
  * Petit run cosmo

Résultat finaux, on a bien formation de structures dans les runs cosmo, mais l'absence de conditions périodiques mène à un effondrement de la boîte et des structures très rapidement.

<img src="https://github.com/mdelorme/gravity_solver_stream/blob/main/figs/cosmo64.png" width="300" height="280" />

## Infos générales

Ces streams sont mensuels. La plupart des communications ont lieu sur mon compte twitter (@astro_md). N'hésitez pas à me contacter pour plus d'informations.
