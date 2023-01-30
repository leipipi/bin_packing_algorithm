
# compile micro kernel
g++ -O2 -c ConfigReader.cpp -o ConfigReader.o
g++ -O2 -fPIC -c Machine.cpp -o Machine.o
g++ -O2 -fPIC -c Timer.cpp -o Timer.o
if [[ $(grep Microsoft /proc/version) ]]; then
        g++ -O2 -c -DWITHOUT_NUMPY -I/usr/include/python2.7 VSSimKernel.cpp -o VSSimKernel.o
        # link micro kernel
        g++ -O2 -DWITHOUT_NUMPY VSSimKernel.o ConfigReader.o Machine.o Timer.o -o VSSimKernel -ldl
else 
        g++ -O0 -g -std=c++11 -c -DWITHOUT_NUMPY -DPLOT -DSAVEJOBS -I/usr/include/python2.7 VSSimKernel.cpp -o VSSimKernel.o
        # link micro kernel
        g++ -O0 -g -std=c++11 -DWITHOUT_NUMPY VSSimKernel.o ConfigReader.o Machine.o Timer.o -o VSSimKernel -lpython2.7 -ldl
fi

# build libraries
g++ -shared -fPIC testlib.cpp -o testlib.so
g++ -shared -fPIC GWFLoader.cpp -o GWFLoader.so
g++ -shared -fPIC SWFLoader.cpp -o SWFLoader.so
g++ -shared -fPIC C0Loader.cpp -o C0Loader.so
g++ -shared -fPIC CXLoader.cpp -o CXLoader.so
g++ -shared -fPIC C6Loader.cpp -o C6Loader.so
g++ -shared -fPIC TestLoader.cpp -o TestLoader.so
g++ -c -fPIC CNS.cpp -o CNS.o
g++ -shared -fPIC CNS.o Machine.o -o CNS.so
g++ -c -fPIC LS.cpp -o LS.o
g++ -shared -fPIC LS.o Machine.o -o LS.so
g++ -c -fPIC LSR.cpp -o LSR.o
g++ -shared -fPIC LSR.o Machine.o -o LSR.so
g++ -c -fPIC SPVS.cpp -o SPVS.o
g++ -shared -fPIC SPVS.o Machine.o -o SPVS.so
g++ -c -fPIC SE.cpp -o SE.o
g++ -shared -fPIC SE.o Machine.o -o SE.so
g++ -c -fPIC KouVS_FF.cpp -o KouVS_FF.o
g++ -shared -fPIC KouVS_FF.o Machine.o -o KouVS_FF.so
g++ -c -fPIC KouVS_FFD_LInf.cpp -o KouVS_FFD_LInf.o
g++ -shared -fPIC KouVS_FFD_LInf.o Machine.o -o KouVS_FFD_LInf.so
g++ -c -fPIC PaniDP.cpp -o PaniDP.o
g++ -shared -fPIC PaniDP.o Machine.o -o PaniDP.so
g++ -c -fPIC PaniL1.cpp -o PaniL1.o
g++ -shared -fPIC PaniL1.o Machine.o -o PaniL1.so
g++ -c -fPIC PaniL2.cpp -o PaniL2.o
g++ -shared -fPIC PaniL2.o Machine.o -o PaniL2.so
g++ -c -fPIC PaniLinf.cpp -o PaniLinf.o
g++ -shared -fPIC PaniLinf.o Machine.o -o PaniLinf.so
g++ -c -fPIC HybridDP.cpp -o HybridDP.o
g++ -shared -fPIC HybridDP.o Machine.o -o HybridDP.so
g++ -c -fPIC HybridL1.cpp -o HybridL1.o
g++ -shared -fPIC HybridL1.o Machine.o -o HybridL1.so
g++ -c -fPIC HybridL2.cpp -o HybridL2.o
g++ -shared -fPIC HybridL2.o Machine.o -o HybridL2.so
g++ -c -fPIC HybridLinf.cpp -o HybridLinf.o
g++ -shared -fPIC HybridLinf.o Machine.o -o HybridLinf.so
g++ -c -fPIC KouVS_FFD_Sum.cpp -o KouVS_FFD_Sum.o
g++ -shared -fPIC KouVS_FFD_Sum.o Machine.o -o KouVS_FFD_Sum.so
g++ -c -fPIC KouVS_FFD_Lexi.cpp -o KouVS_FFD_Lexi.o
g++ -shared -fPIC KouVS_FFD_Lexi.o Machine.o -o KouVS_FFD_Lexi.so
g++ -c -fPIC KouVS_FFD_Lexi_Reordered.cpp -o KouVS_FFD_Lexi_Reordered.o
g++ -shared -fPIC KouVS_FFD_Lexi_Reordered.o Machine.o -o KouVS_FFD_Lexi_Reordered.so
g++ -c -fPIC KouVS_FFD_L2.cpp -o KouVS_FFD_L2.o
g++ -shared -fPIC KouVS_FFD_L2.o Machine.o -o KouVS_FFD_L2.so
g++ -c -fPIC KouVS_BF.cpp -o KouVS_BF.o
g++ -shared -fPIC KouVS_BF.o Machine.o -o KouVS_BF.so
g++ -c -fPIC KouVS_BFD_LInf.cpp -o KouVS_BFD_LInf.o
g++ -shared -fPIC KouVS_BFD_LInf.o Machine.o -o KouVS_BFD_LInf.so
g++ -c -fPIC KouVS_BFD_Sum.cpp -o KouVS_BFD_Sum.o
g++ -shared -fPIC KouVS_BFD_Sum.o Machine.o -o KouVS_BFD_Sum.so
g++ -c -fPIC KouVS_BFD_Lexi.cpp -o KouVS_BFD_Lexi.o
g++ -shared -fPIC KouVS_BFD_Lexi.o Machine.o -o KouVS_BFD_Lexi.so
g++ -c -fPIC KouVS_BFD_Lexi_Reordered.cpp -o KouVS_BFD_Lexi_Reordered.o
g++ -shared -fPIC KouVS_BFD_Lexi_Reordered.o Machine.o -o KouVS_BFD_Lexi_Reordered.so
g++ -c -fPIC KouVS_BFD_L2.cpp -o KouVS_BFD_L2.o
g++ -shared -fPIC KouVS_BFD_L2.o Machine.o -o KouVS_BFD_L2.so
