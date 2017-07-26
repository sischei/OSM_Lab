rm noVec
rm Vec
echo  "======================================================================="
echo  "Running noVec ..."
g++ -O2  main.cpp -o noVec
./noVec
echo  "======================================================================="
echo  "Running Vec ...."
g++ -O2 -ftree-vectorize -ftree-vectorizer-verbose=2 -fopt-info-vec main.cpp -o Vec
./Vec
