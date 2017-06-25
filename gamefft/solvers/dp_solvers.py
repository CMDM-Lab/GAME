"""
File name  : dp_solvers.py 
Description:
Author     : Alioune Schurz
Lab        : Computational Molecular Design and Detection Lab (CMDD lab)
Institute  : National Taiwan University
"""
from ctypes import *
import math
from numpy import *
import gamefft.npsdb.data_handler as dat
import gamefft.solvers.template as template
import gamefft.tools 
import os.path

# Abstract classes ===============================
class AbstractDPSolver(template.AbstractSolver):
        
    def __init__(self):
        template.AbstractSolver.__init__(self)
    # @abstract
    def run(self, csccp, run_options={}): 
        return dat.CSCCP.Result(None) 
    
    def backtrack_selected_sidechain(self, s, w, r, W, L):
        """
        Backtracks the sidechain info list to obtain the corresponding compound 
        of the given s,w,r 
        Args:
            s    : chosen position index
            w    : chosen target weight
            W    : Weight list (2D array)
            L    : Selected sidechains info list
        returns: 
            found sidechain (1D array)
        """    
        if s == -1:
            return []
        else:
            j = L[s][w][r][0]
            t = L[s][w][r][1]
            return self.backtrack_selected_sidechain(s - 1, w - W[s][j], t, W, L) + [j]
    
    def compute_real_weight(self, csccp, w, r, W, L):
        comp = self.backtrack_selected_sidechain(csccp.n - 1, w, r, W, L)
        weights = [0] + [csccp.W[s][comp[s]] for s in range(0, csccp.n)]
        return sum(weights), comp

class AbstractDPSolverCompressed(template.AbstractSolver):
        
    def __init__(self):
        self.float_t = float64
        self.int_t = int64
        
    # @abstract
    def run(self, csccp, run_options={}): 
        return dat.CSCCP.Result(None) 
    
    def backtrack_selected_sidechain(self, s, w, r, W, L, K_encoder, R_encoder):
        """
        Backtracks the sidechain info list to obtain the corresponding compound 
        of the given s,w,r 
        Args:
            s    : chosen position index
            w    : chosen target weight
            W    : Weight list (2D array)
            L    : Selected sidechains info list
        returns: 
            found sidechain (1D array)
        """    
        if s == -1:
            return []
        else:
            j = K_encoder.decode_at(L[w][r][0], s)  
            t = R_encoder.decode_at(L[w][r][1], s)
            return self.backtrack_selected_sidechain(s - 1, w - W[s][j], t, W, L, K_encoder, R_encoder) + [j]
    
    def compute_real_weight(self, csccp, w, r, W, L, K_encoder, R_encoder):
        comp = self.backtrack_selected_sidechain(csccp.n - 1, w, r, W, L, K_encoder, R_encoder)
        weights = [0] + [csccp.W[s][comp[s]] for s in range(0, csccp.n)]
        return sum(weights), comp

class AbstractDPKernel:
    # @abstract
    def run(self, n, R, wt_max, K, W, P, L, C, run_options):
        SIGINT, TIMEOUT = (0, -1)
        return (SIGINT, TIMEOUT)

class AbstractDPKernelCompressed:
    # @abstract
    def run(self, n, R, wt_max, K, W, P, L, C, K_encoder, R_encoder, run_options):
        SIGINT, TIMEOUT = (0, -1)
        return (SIGINT, TIMEOUT)

class AbstractEncoder:
    # @abstract
    def encode(self):
        return
    # @abstract
    def decode(self, csccp, w, r, W, L):
        return
        
# Encoders =======================================

class Encoder1(AbstractEncoder):
    
    class NotEncodableException(Exception):
        def __init__(self, message):
            Exception.__init__(self, message)

    def __init__(self):
        self.K = None
        self.Km1 = None
        self.Kbitlen = None
        self.offsets = None
                
    def setup(self, K, maxbits):
        self.K = K
        self.Km1 = [k - 1 for k in self.K]
        if not self.encodable(self.Km1, maxbits):
            raise self.NotEncodableException('Not encodable in %d bits.' % maxbits)
        self.Kbitlen = [self.__bitlen(k) for k in self.Km1]
        self.offsets = self.__get_offsets(self.Kbitlen)
        
    def encodable(self, K, maxbits):
        s = sum([self.__bitlen(k) for k in K])
        return s <= maxbits
    
    def encode_at(self, result, val, i):
        return self.__encode_item(self.offsets, result, val, i)
    
    def decode_at(self, result, i):
        return self.__decode_item(self.offsets, self.Kbitlen, result, i)
    
    def encode_all(self, encoded):
        return self.__encode_all(self.offsets, 0, encoded)
    
    def decode_all(self, result):
        return self.__decode_all(self.offsets, self.Km1, self.Kbitlen, result)
    
    def __bitlen(self, n):
        if n==0:
            return 0
        else:
            return int(math.floor(math.log(n, 2)) + 1)
    
    def __write_at(self, a, b, offset):
        return uint32(a | b << offset)
    
    def __get_at(self, a, offset, width):
        return uint32((a & (2 ** width - 1) << offset) >> offset)

    def __cumulative_sum(self, Val):
        ValAcc = [0 for x in Val]
        ValAcc[0] = Val[0]
        for i in range(1, len(Val)):
            ValAcc[i] = Val[i] + ValAcc[i - 1]
        return ValAcc

    def __get_offsets(self, Kbitlen):
        offsets = [0] + self.__cumulative_sum(Kbitlen[:-1])   
        return offsets
                
    def __encode_item(self, offsets, result, val, i):
        return self.__write_at(result, val, offsets[i])
    
    def __decode_item(self, offsets, Kbitlen, result, i):
        return self.__get_at(result, offsets[i], Kbitlen[i])
    
    def __encode_all(self, offsets, result, val):
        for i, v in enumerate(val):
            result = self.__encode_item(offsets, result, v, i)
        return result
            
    def __decode_all(self, offsets, Km1, Kbitlen, result):
        res = []
        for i, k in enumerate(Km1):
            res.append(self.__decode_item(offsets, Kbitlen, result, i))
        return res  
    
    # ==========Don't touch========
    def decode(self, s, w, r, W, L):
        if s == -1:
            return []
        else:
            j = L[s][w][r][0]
            t = L[s][w][r][1]
            return self.decode(s - 1, w - W[s][j], t, W, L) + [j]
    
    def compute_real_weight(self, csccp, w, r, K, W, L): 
        comp = self.decode(csccp.n - 1, w, r, W, L)
        weights = [0] + [csccp.W[s][comp[s]] for s in range(0, csccp.n)]
        return sum(weights), comp
    
# Kernels ========================================

class DPKernelCpuPy(AbstractDPKernel):
    def run(self, n, R, wt_max, K, W, P, L, C, run_options):
        for s in range(0, n):
            for w in range(1, wt_max + 1):
                A = []
                for k in range(0, K[s]):
                    if w - W[s][k] < 0:
                        continue
                    for r in range(0, R):
                        if s > 0:
                            A.append((k, r, P[s][k] * C[s - 1][w - W[s][k]][r]))
                        else:
                            if r == 1 and w - W[s][k] == 0:
                                A.append((k, r, P[s][k]))
                    # -- END for r 
                # -- END for k
                # sort by probability
                A = sorted(A, key=lambda x:x[2], reverse=True)
                # update C and L
                for r in range(0, min(R, len(A))):
                    L[s][w][r][0] = A[r][0]
                    L[s][w][r][1] = A[r][1]
                    C[s][w][r] = A[r][2]
            # -- END for w
        # -- END for s
        SIGINT, TIMEOUT = 0, -1
        return (SIGINT, TIMEOUT)

class DPKernelCpuPyOpt(AbstractDPKernel):         
    def run(self, n, R, wt_max, K, W, P, L, C, run_options):
        Cnew = copy(C)
        Cswi = [C, Cnew]
        for s in range(0, n):
            prevIdx = s % 2
            newIdx = (s + 1) % 2
            for w in range(1, wt_max + 1):
                A = []
                for k in range(0, K[s]):
                    if w - W[s][k] < 0:
                        continue
                    for r in range(0, R):
                        if s > 0:
                            A.append((k, r, P[s][k] * Cswi[prevIdx][w - W[s][k]][r]))
                        else:
                            if r == 1 and w - W[s][k] == 0:
                                A.append((k, r, P[s][k]))
                    # -- END for r 
                # -- END for k
                # sort by probability
                A = sorted(A, key=lambda x:x[2], reverse=True)
                # update C and L
                for r in range(0, min(R, len(A))):
                    L[s][w][r][0] = A[r][0]
                    L[s][w][r][1] = A[r][1]
                    Cswi[newIdx][w][r] = A[r][2]
            # -- END for w   
        # -- END for s
        C = Cswi[(s + 1) % 2]
        SIGINT, TIMEOUT = 0, -1
        return (SIGINT, TIMEOUT)

class DPKernelCpuCcOpt(AbstractDPKernel):
    def run(self, n, R, wt_max, K, W, P, L, C, run_options):
        dll_name = "libgamefft.so"
        myDll = cdll.LoadLibrary(dll_name)
        kernelCcOpt = myDll.DPKernelCpuCcOpt
        kernelCcOpt.argtypes = [c_int, c_int, c_int,
                                c_int, ctypeslib.ndpointer(c_int), ctypeslib.ndpointer(c_int), ctypeslib.ndpointer(c_float),
                                ctypeslib.ndpointer(c_uint16), ctypeslib.ndpointer(c_float)]
        l = max([len(x) for x in W])
        TIMEOUT = kernelCcOpt(n, R, wt_max, l, K, W, P, L, C)  
        SIGINT = 0
        return (SIGINT, TIMEOUT)   
    
class DPKernelCpuCcOptBubble(AbstractDPKernel):
    def run(self, n, R, wt_max, K, W, P, L, C, run_options):
        dll_name = "libgamefft.so"
        myDll = cdll.LoadLibrary(dll_name)
        kernelCcOpt = myDll.DPKernelCpuCcOptBubble
        kernelCcOpt.argtypes = [c_int, c_int, c_int,
                                c_int, ctypeslib.ndpointer(c_int), ctypeslib.ndpointer(c_int), ctypeslib.ndpointer(c_float),
                                ctypeslib.ndpointer(c_uint16), ctypeslib.ndpointer(c_float)]
        l = max([len(x) for x in W])
        TIMEOUT = kernelCcOpt(n, R, wt_max, l, K, W, P, L, C)  
        SIGINT = 0
        return (SIGINT, TIMEOUT)  
             
class DPKernelCpuPyOptCompressed(AbstractDPKernelCompressed): 
               
    def run(self, n, R, wt_max, K, W, P, L, C, K_encoder, R_encoder, run_options):

        Cnew = copy(C)
        Cswi = [C, Cnew]

        for s in range(0, n):
            prevIdx = s % 2
            newIdx = (s + 1) % 2
            for w in range(1, wt_max + 1):
                A = []
                for k in range(0, K[s]):
                    if w - W[s][k] < 0:
                        continue
                    for r in range(0, R):
                        if s > 0:
                            A.append((k, r, P[s][k] * Cswi[prevIdx][w - W[s][k]][r]))
                        else:
                            if r == 1 and w - W[s][k] == 0:
                                A.append((k, r, P[s][k]))
                    # -- END for r 
                # -- END for k
                # sort by probability
                A = sorted(A, key=lambda x:x[2], reverse=True)
                # update C and L
                for r in range(0, min(R, len(A))):
                    # L[s][w][r][0] = A[r][0]
                    # L[s][w][r][1] = A[r][1]
                    L[w][r][0] = K_encoder.encode_at(L[w][r][0], A[r][0], s)
                    L[w][r][1] = R_encoder.encode_at(L[w][r][1], A[r][1], s)
                    Cswi[newIdx][w][r] = A[r][2]
            # -- END for w   
        # -- END for s
        C = Cswi[(s + 1) % 2]

        SIGINT, TIMEOUT = 0, -1
        return (SIGINT, TIMEOUT)
    
class DPKernelCpuCcOptBubbleCompressed(AbstractDPKernelCompressed):
    def run(self, n, R, wt_max, K, W, P, L, C, K_encoder, R_encoder, run_options):
        K_offsets = array(K_encoder.offsets, dtype=int32)
        R_offsets = array(R_encoder.offsets, dtype=int32)
        myDll = cdll.LoadLibrary("libgamefft.so")
        kernelCcOptCompressed = myDll.DPKernelCpuCcOptBubbleCompressed
        kernelCcOptCompressed.argtypes = [c_int, c_int, c_int,
                                c_int, ctypeslib.ndpointer(c_int), ctypeslib.ndpointer(c_int), ctypeslib.ndpointer(c_float),
                                ctypeslib.ndpointer(c_uint32), ctypeslib.ndpointer(c_float),
                                ctypeslib.ndpointer(c_int), ctypeslib.ndpointer(c_int)]
                
        l = max([len(x) for x in W])
        TIMEOUT = kernelCcOptCompressed(n, R, wt_max, l, K, W, P, L, C, K_offsets, R_offsets)  

        SIGINT = 0
        return (SIGINT, TIMEOUT)          

# Generic solvers =================================

class GenericDPSolver(AbstractDPSolver):
    def __init__(self, kernel):
        AbstractDPSolver.__init__(self)
        self.kernel = kernel

    def run(self, csccp, run_options):   
        if run_options["verbose"]:
            print "Starting Dynamic Programming"
        if csccp.trivial():
            raise dat.TrivialCSCCP("This CSCCP problem is trivial.")
    
        # 1- Variable init   * * * * * * * * * * * * * * * * * * * * * * * *
        l = max([len(x) for x in csccp.W])
        n = csccp.n
        R = csccp.R
        D = run_options["dec"]
        wt_max = int(csccp.wmax * 10 ** D + n * 0.5 / 10 ** D)
        wt_min = int(csccp.wmin * 10 ** D - n * 0.5 / 10 ** D)
        W = array([([int(w * 10 ** D) for w in wside] + [0] * l)[:l] for wside in csccp.W])
        P = array([(pside + [0] * l)[:l] for pside in csccp.P])
        K = array(csccp.K)
        C = zeros((n, wt_max + 1, R), self.float_t)  # [0 .. n-1] x [0 .. wt_max] x [0 .. R-1] 
        L = zeros((n, wt_max + 1, R, 2), self.int_t)
        
        try:
            # 2- Core computations   * * * * * * * * * * * * * * * * * * * *  
            SIGINT, TIMEOUT = self.kernel.run(n, R, wt_max, K, W, P, L, C, run_options)
            if SIGINT:    
                raise KeyboardInterrupt
            if TIMEOUT > 0:
                raise tools.Timer.TimeoutException(TIMEOUT)
            # 3- Result filtering starting * * * * * * * * * * * * * * * * *
            valid_candidates = []
            len_filtered_results = 0
            for w in range(wt_min, wt_max + 1):
                for r in range(0, R):
                    probability = C[n - 1][w][r]
                    if probability > 0:
                        len_filtered_results += 1
                        real_weight, comp = self.compute_real_weight(csccp, w, r, W, L)
                        if real_weight >= csccp.wmin and real_weight <= csccp.wmax:
                            valid_candidates.append((probability * csccp.scaffold_probability_rel2config, real_weight + csccp.scaffold_weight_rel2config, comp))        
    
            valid_candidates = sorted(valid_candidates, key=lambda x:x[0], reverse=True)
            valid_candidates = valid_candidates[0:R]
            return dat.CSCCP.Result(array(valid_candidates, dtype=dtype([('p', float), ('w', float), ('comp', list)])))
        except KeyboardInterrupt:
            print "*  WARNING: the program was interrupted by user (CTRL+C). The program will now terminate *"
            exit()

class GenericDPSolverOpt(AbstractDPSolver):
    
    def __init__(self, kernelOpt):
        AbstractDPSolver.__init__(self)
        self.kernelOpt = kernelOpt

    def run(self, csccp, run_options):
                
        if run_options["verbose"]:
            print "Starting Dynamic Programming"
        if csccp.trivial():
            raise dat.TrivialCSCCP("This CSCCP problem is trivial.")
    
        # 1- Variable init   * * * * * * * * * * * * * * * * * * * * * * * *
        l = max([len(x) for x in csccp.W])
        n = csccp.n
        R = csccp.R
        D = run_options["dec"]
        wt_max = int(csccp.wmax * 10 ** D + n * 0.5 / 10 ** D)
        wt_min = int(csccp.wmin * 10 ** D - n * 0.5 / 10 ** D)
        W = array([([int(w * 10 ** D) for w in wside] + [0] * l)[:l] for wside in csccp.W], dtype=int32)
        P = array([(pside + [0] * l)[:l] for pside in csccp.P], dtype=float32)
        K = array(csccp.K, dtype=int32)
        C = zeros((wt_max + 1, R), self.float_t)  # [0 .. n-1] x [0 .. wt_max] x [0 .. R-1] 
        L = zeros((n, wt_max + 1, R, 2), self.int_t)
        
        try:
            # 2- Core computations   * * * * * * * * * * * * * * * * * * * *  
            SIGINT, TIMEOUT = self.kernelOpt.run(n, R, wt_max, K, W, P, L, C, run_options)
            if SIGINT:    
                raise KeyboardInterrupt
            if TIMEOUT > 0:
                raise tools.Timer.TimeoutException(TIMEOUT)
            # 3- Result filtering starting * * * * * * * * * * * * * * * * *
            valid_candidates = []
            len_filtered_results = 0
            for w in range(wt_min, wt_max + 1):
                for r in range(0, R):
                    probability = C[w][r]
                    if probability > 0:
                        len_filtered_results += 1
                        real_weight, comp = self.compute_real_weight(csccp, w, r, W, L)
                        if real_weight >= csccp.wmin and real_weight <= csccp.wmax:
                            valid_candidates.append((probability * csccp.scaffold_probability_rel2config, real_weight + csccp.scaffold_weight_rel2config, comp))        
    
            valid_candidates = sorted(valid_candidates, key=lambda x:x[0], reverse=True)
            valid_candidates = valid_candidates[0:R]
            return dat.CSCCP.Result(array(valid_candidates, dtype=dtype([('p', float), ('w', float), ('comp', list)])))
        except KeyboardInterrupt:
            print "*  WARNING: the program was interrupted by user (CTRL+C). The program will now terminate *"
            exit()

class GenericDPSolverOptCompressed(AbstractDPSolverCompressed):
    
    def __init__(self, kernelOptCompressed):
        AbstractDPSolverCompressed.__init__(self)
        self.kernelOptCompressed = kernelOptCompressed
        self.int_t = uint32
        self.float_t = float32

    def run(self, csccp, run_options):  
        if run_options["verbose"]:
            print "Starting Dynamic Programming C-optimized L-compressed"
        if csccp.trivial():
            raise dat.TrivialCSCCP("This CSCCP problem is trivial.")
    
        # 1- Variable init   * * * * * * * * * * * * * * * * * * * * * * * *
        l = max([len(x) for x in csccp.W])
        n = csccp.n
        R = csccp.R
        D = run_options["dec"]
        wt_max = int(csccp.wmax * 10 ** D + n * 0.5 / 10 ** D)
        wt_min = int(csccp.wmin * 10 ** D - n * 0.5 / 10 ** D)
        W = array([([int(w * 10 ** D) for w in wside] + [0] * l)[:l] for wside in csccp.W], dtype=int32)
        P = array([(pside + [0] * l)[:l] for pside in csccp.P], dtype=float32)
        K = array(csccp.K, dtype=int32)
        C = zeros((wt_max + 1, R), float32)  # [0 .. n-1] x [0 .. wt_max] x [0 .. R-1] 
        L = zeros((wt_max + 1, R, 2), uint32)
        K_encoder = Encoder1()
        R_encoder = Encoder1()
        K_encoder.setup(K, 32)        
        R_encoder.setup([R - 1 for x in range(0, n)], 32)
        
        try:
            # 2- Core computations   * * * * * * * * * * * * * * * * * * * *  
            SIGINT, TIMEOUT = self.kernelOptCompressed.run(n, R, wt_max, K, W, P, L, C, K_encoder, R_encoder, run_options)
                    
            if SIGINT:    
                raise KeyboardInterrupt
            if TIMEOUT > 0:
                raise tools.Timer.TimeoutException(TIMEOUT)
            # 3- Result filtering starting * * * * * * * * * * * * * * * * *
            valid_candidates = []
            len_filtered_results = 0
            for w in range(wt_min, wt_max + 1):
                for r in range(0, R):
                    probability = C[w][r]
                    if probability > 0:
                        len_filtered_results += 1
                        real_weight, comp = self.compute_real_weight(csccp, w, r, W, L, K_encoder, R_encoder)
                        if real_weight >= csccp.wmin and real_weight <= csccp.wmax:
                            valid_candidates.append((probability * csccp.scaffold_probability_rel2config, real_weight + csccp.scaffold_weight_rel2config, comp))        
    
            valid_candidates = sorted(valid_candidates, key=lambda x:x[0], reverse=True)
            valid_candidates = valid_candidates[0:R]
            return dat.CSCCP.Result(array(valid_candidates, dtype=dtype([('p', float), ('w', float), ('comp', list)])))
        except KeyboardInterrupt:
            print "*  WARNING: the program was interrupted by user (CTRL+C). The program will now terminate *"
            exit()
        
# Specialized DP solvers =============================

class DPSolverCpuPy(GenericDPSolver):
    def __init__(self):
        GenericDPSolver.__init__(self, DPKernelCpuPy())

class DPSolverCpuPyOpt(GenericDPSolverOpt):
    def __init__(self):
        GenericDPSolverOpt.__init__(self, DPKernelCpuPyOpt())

class DPSolverCpuCcOpt(GenericDPSolverOpt):
    def __init__(self):
        GenericDPSolverOpt.__init__(self, DPKernelCpuCcOpt())
        
class DPSolverCpuPyOptCompressed(GenericDPSolverOptCompressed):
    def __init__(self):
        GenericDPSolverOptCompressed.__init__(self, DPKernelCpuPyOptCompressed())
    
class DPSolverCpuCcOptBubbleCompressed(GenericDPSolverOptCompressed):
    def __init__(self):
        GenericDPSolverOptCompressed.__init__(self, DPKernelCpuCcOptBubbleCompressed())        
    
