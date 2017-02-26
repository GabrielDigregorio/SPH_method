#Shell script for performance testing
echo "N influence for equally spaced particules"
./Neighboor_performance.exe 1 1.2 3 0
./Neighboor_performance.exe 1 1.2 5 0
./Neighboor_performance.exe 1 1.2 7 0
./Neighboor_performance.exe 1 1.2 10 0
./Neighboor_performance.exe 1 1.2 12 0
./Neighboor_performance.exe 1 1.2 15 0
./Neighboor_performance.exe 1 1.2 18 0


echo "kh influence for equally spaced particules, N = 14^3"
./Neighboor_performance.exe 1 1.2 13 0
./Neighboor_performance.exe 1 1.5 13 0
./Neighboor_performance.exe 1 1.7 13 0
./Neighboor_performance.exe 1 2.3 13 0
./Neighboor_performance.exe 1 2.5 13 0
./Neighboor_performance.exe 1 2.7 13 0
./Neighboor_performance.exe 1 3.0 13 0

echo "N influence for equally spaced particules"
./Neighboor_performance.exe 1 1.2 3 0.2
./Neighboor_performance.exe 1 1.2 5 0.2
./Neighboor_performance.exe 1 1.2 7 0.2
./Neighboor_performance.exe 1 1.2 10 0.2
./Neighboor_performance.exe 1 1.2 12 0.2
./Neighboor_performance.exe 1 1.2 15 0.2
./Neighboor_performance.exe 1 1.2 18 0.2


echo "kh influence for equally spaced particules, N = 14^3"
./Neighboor_performance.exe 1 1.2 13 0.2
./Neighboor_performance.exe 1 1.5 13 0.2
./Neighboor_performance.exe 1 1.7 13 0.2
./Neighboor_performance.exe 1 2.3 13 0.2
./Neighboor_performance.exe 1 2.5 13 0.2
./Neighboor_performance.exe 1 2.7 13 0.2
./Neighboor_performance.exe 1 3.0 13 0.2
