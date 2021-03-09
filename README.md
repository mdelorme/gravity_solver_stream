# Gravitational solver stream

Dépôt github correspondant aux streams de développements de ma chaîne twitch (https://www.twitch.tv/astro_md).

## Session 1 (11/12/2020)
Lien vidéo : https://www.youtube.com/watch?v=dXc-2P8fHfk

Contenu : 
 * Explication du problème à n-corps. 
 * Présentation de la méthode de sommation directe.
 * Implémentation d'un système à 2 corps, Terre-Soleil
 <img src="https://github.com/mdelorme/gravity_solver_stream/blob/main/figs/session1.png" width="300" height="300" />
 
## Session 2 (15/01/2021) 
Lien vidéo : https://www.youtube.com/watch?v=__KIreszA6I
 
Contenu :
 * Nettoyage du code pour une approche objet
 * Extension à un nombre arbitraire de particules
 * Mise en place de particules traceuses 
 * Exemple jouet sur la collision de deux galaxies
 
 <img src="https://github.com/mdelorme/gravity_solver_stream/blob/main/figs/session2.gif" width="300" height="300" />
 
## Session 3 (05/03/2021)
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

## Session 4 (Date à décider)

Contenu :
  * Optimisation du code pour GPU
  * Introduction de parallélisme hiérarchique
  * Run plus avec un plus grand nombre de particules massives
 
## Infos générales

Ces streams sont mensuels. La plupart des communications ont lieu sur mon compte twitter (@astro_md). N'hésitez pas à me contacter pour plus d'informations.
