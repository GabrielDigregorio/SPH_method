#Shell script for performance testing
echo "N influence for equally spaced particules"
./sph 1 1.2 3 0
./sph 1 1.2 5 0
./sph 1 1.2 7 0
./sph 1 1.2 10 0
./sph 1 1.2 12 0
./sph 1 1.2 15 0
./sph 1 1.2 18 0


echo "kh influence for equally spaced particules, N = 14^3"
./sph 1 1.2 13 0
./sph 1 1.5 13 0
./sph 1 1.7 13 0
./sph 1 2.3 13 0
./sph 1 2.5 13 0
./sph 1 2.7 13 0
./sph 1 3.0 13 0

echo "N influence for equally spaced particules"
./sph 1 1.2 3 0.2
./sph 1 1.2 5 0.2
./sph 1 1.2 7 0.2
./sph 1 1.2 10 0.2
./sph 1 1.2 12 0.2
./sph 1 1.2 15 0.2
./sph 1 1.2 18 0.2


echo "kh influence for equally spaced particules, N = 14^3"
./sph 1 1.2 13 0.2
./sph 1 1.5 13 0.2
./sph 1 1.7 13 0.2
./sph 1 2.3 13 0.2
./sph 1 2.5 13 0.2
./sph 1 2.7 13 0.2
./sph 1 3.0 13 0.2
