import numpy as np

class Ephem_Struct:
        
    def __init__(self):

        self.EPHEM_ARRAY_SIZE = 9 # сколько один и тот же спутник будет встречаться в RINEX файле
        self.PRN_IGSO = [6, 7, 8, 9, 10, 13, 16, 38, 39, 40] # номера спутников на IGSO

        # Initialize the Ephemeris data before decoding.

        self.PRN    = None
        self.SatH1  = None
        self.AODC   = None
        self.URAI   = None
        self.WN     = None

        self.T_GD_1 = None
        self.T_GD_2 = None

        # Ionospheric Delay Model Parameters (alpha, beta)

        self.alpha0 = None
        self.alpha1 = None
        self.alpha2 = None
        self.alpha3 = None

        self.beta0  = None
        self.beta1  = None
        self.beta2  = None
        self.beta3  = None

        # Clock Correction Parameters

        self.a0     = None
        self.a1     = None
        self.a2     = None
        self.t_oc   = None

        self.a0_utc = None
        self.a1_utc = None
        self.t_ref  = None
        self.w_ref  = None
        self.dtls   = None

        # Age of Data, Ephemeris (IODE)

        self.AODE   = None

        # Ephemeris Parameters

        self.deltan = None
        self.C_uc   = None
        self.M_0    = None
        self.e      = None

        self.C_us   = None
        self.C_rc   = None
        self.C_rs   = None

        self.sqrtA    = None
        self.i_0      = None
        self.C_ic     = None
        self.omegaDot = None
        self.C_is     = None
        self.iDot     = None
        self.omega_0  = None
        self.omega    = None
        self.t_oe     = None

        # Tow of first decoded 
            
        self.SOW      = None

        # Time parametrs

        self.y        = None
        self.m        = None
        self.d        = None
        self.hh       = None
        self.mm       = None
        self.sec      = None

    def rename_exp_designator(self, x):
        return x.replace('D', 'E')
    
    def date2bdt(self, year, mounth, day, hour, minute, second):

        day_in_year = [0,31,59,90,120,151,181,212,243,273,304,334]

        dyear = year - 2006 # сколько прошло лет с 2006

        lpdays = int(dyear/4) + 1 # сколько високосных дней нужно добавить
        if (dyear%2 == 0) and (mounth <= 2):
            lpdays -= 1

        dn = dyear * 365 + day_in_year[mounth - 1] + day + lpdays - 1
        wn = int(dn/7)
        sec = (dn%7) * 86400 + hour * 3600 + minute * 60 + second

        return dn, wn, sec 
    
    def subBdsTime(self, ref_wn, tmp_wn, ref_time, tmp_time, dt = 0):
        
        dt = tmp_time - ref_time
        dt += (tmp_wn - ref_wn) * 604800

        return dt
    
    def read_from_lastSTR(self, file, str_lines, num_str, max_char):

        str_lines = str_lines[(max_char + 1)*(num_str-1):]
        str_lines = file.read(max_char+1)

        return str_lines[1:]

    def pars_rinex(self, MAX_SAT, file_name, max_char, num_rinex_str, mode = 1):

        eph = np.array([[]], dtype='float') #
        all_eph = np.array([[]], dtype='float') #
        
        NUM_SAT = 0 # количество спутников, которое поместилось в массив размером (ieph, EPHEM_ARRAY_SIZE)
        health_flag = 0 # флаг указывающий на то, что эфемериды для данного спутника прочитана

        file = open(str(file_name) + '.txt', 'r')

        # Read Header 

        while(1):
                    
            rinex_str = file.read(max_char)

            if (rinex_str == None):
                break
                          
            if (rinex_str[60:] == 'END OF HEADER'):
                break
                    
            elif (rinex_str[:5] == 'BDSA'):
                self.alpha0 = float(self.rename_exp_designator(rinex_str[6:17]))
                self.alpha1 = float(self.rename_exp_designator(rinex_str[18:29]))
                self.alpha2 = float(self.rename_exp_designator(rinex_str[30:41]))
                self.alpha3 = float(self.rename_exp_designator(rinex_str[42:53]))

            elif (rinex_str[:5] == 'BDSB'):
                self.beta0 = float(self.rename_exp_designator(rinex_str[9:17]))
                self.beta1 = float(self.rename_exp_designator(rinex_str[18:29]))
                self.beta2 = float(self.rename_exp_designator(rinex_str[30:41]))
                self.beta3 = float(self.rename_exp_designator(rinex_str[42:53]))

            elif (rinex_str[:5] == 'BDUT'):
                self.a0_utc = float(self.rename_exp_designator(rinex_str[5:22]))
                self.a1_utc = float(self.rename_exp_designator(rinex_str[22:38]))
                self.t_ref = int(rinex_str[43:45])
                self.w_ref = int(rinex_str[47:50])

            elif (rinex_str[:5] == 'LEAP SECONDS'):
                self.dtls = int(rinex_str[5:6])

        # Read Ephemeris

        while(1):

            all_eph = np.append(all_eph, eph, axis=mode)

            if (rinex_str == None):
                print('Пустой файл')
                break

            if health_flag == 1: # сброс флага "здоровья"
                health_flag = 0

            if NUM_SAT == MAX_SAT: # выход из цикла при заполнении массива
                break
            
            rinex_str = file.read(max_char + 81*(num_rinex_str - 1)) # чтение строки, следующей за последней прочитанной
            if num_rinex_str > 1:
                rinex_str = self.read_from_lastSTR(file, rinex_str, num_rinex_str, max_char)

            if rinex_str[:1] != 'C':
                print('Файл эфемерид не соответствует BDS')
                break

            self.PRN = int(rinex_str[1:3])

            if self.PRN in self.PRN_IGSO:

                self.y = int(rinex_str[4:8])
                self.m = int(rinex_str[9:11])
                self.d = int(rinex_str[12:14])

                self.hh = int(rinex_str[15:17])
                self.mm = int(rinex_str[18:20])
                self.sec = int(rinex_str[21:23])
                
                num_week, num_sec = self.date2bdt(self.y, self.m, self.d, self.hh, self.mm, self.sec)
                
                # Массив с количеством строк равным количеству спутников на IGSO

                # EPOCH, SV CLK

                eph = np.append(eph, num_week, num_sec) 

                self.a0 = float(self.rename_exp_designator(rinex_str[23:42]))
                eph = np.append(eph, self.a0) 

                self.a1 = float(self.rename_exp_designator(rinex_str[42:61]))
                eph = np.append(eph, self.a1) 

                self.a2 = float(self.rename_exp_designator(rinex_str[61:80]))
                eph = np.append(eph, self.a2) 

                # BROADCAST ORBIT-1

                rinex_str = file.read(max_char)
                num_rinex_str += 1

                if (rinex_str == None):
                    break
                
                self.AODE = float(self.rename_exp_designator(rinex_str[4:23]))
                eph = np.append(eph, self.AODE) 

                self.C_rs = float(self.rename_exp_designator(rinex_str[23:42]))
                eph = np.append(eph, self.C_rs) 

                self.deltan = float(self.rename_exp_designator(rinex_str[42:61]))
                eph = np.append(eph, self.deltan)

                self.M_0 = float(self.rename_exp_designator(rinex_str[61:80]))
                eph = np.append(eph, self.M_0)

                # BROADCAST ORBIT-2

                rinex_str = file.read(max_char)
                num_rinex_str += 1

                if (rinex_str == None):
                    break

                self.C_uc = float(self.rename_exp_designator(rinex_str[4:23]))
                eph = np.append(eph, self.C_uc) 

                self.e = float(self.rename_exp_designator(rinex_str[23:42]))
                eph = np.append(eph, self.e) 

                self.C_us = float(self.rename_exp_designator(rinex_str[42:61]))
                eph = np.append(eph, self.C_us)

                self.sqrtA = float(self.rename_exp_designator(rinex_str[61:80]))
                eph = np.append(eph, self.sqrtA)

                # BROADCAST ORBIT-3

                rinex_str = file.read(max_char)
                num_rinex_str += 1

                if (rinex_str == None):
                    break
                
                self.t_oe = float(self.rename_exp_designator(rinex_str[4:23]))
                eph = np.append(eph, self.t_oe) 

                self.C_ic = float(self.rename_exp_designator(rinex_str[23:42]))
                eph = np.append(eph, self.C_ic) 

                self.omega_0 = float(self.rename_exp_designator(rinex_str[42:61]))
                eph = np.append(eph, self.omega_0)

                self.C_is = float(self.rename_exp_designator(rinex_str[61:80]))
                eph = np.append(eph, self.C_is)

                # BROADCAST ORBIT-4

                rinex_str = file.read(max_char)
                num_rinex_str += 1

                if (rinex_str == None):
                    break
                
                self.i_0 = float(self.rename_exp_designator(rinex_str[4:23]))
                eph = np.append(eph, self.i_0) 

                self.C_rc = float(self.rename_exp_designator(rinex_str[23:42]))
                eph = np.append(eph, self.C_rc) 

                self.omega = float(self.rename_exp_designator(rinex_str[42:61]))
                eph = np.append(eph, self.omega)

                self.omegaDot = float(self.rename_exp_designator(rinex_str[61:80]))
                eph = np.append(eph, self.omegaDot)

                # BROADCAST ORBIT-5

                rinex_str = file.read(max_char)
                num_rinex_str += 1

                if (rinex_str == None):
                    break
                
                self.iDot = float(self.rename_exp_designator(rinex_str[4:23]))
                eph = np.append(eph, self.iDot) 

                self.WN = float(self.rename_exp_designator(rinex_str[42:61]))
                eph = np.append(eph, self.WN) 

                # BROADCAST ORBIT-6

                rinex_str = file.read(max_char)
                num_rinex_str += 1

                if (rinex_str == None):
                    break

                self.URAI = float(self.rename_exp_designator(rinex_str[4:23]))
                eph = np.append(eph, self.URAI) 

                self.SatH1 = float(self.rename_exp_designator(rinex_str[23:42]))
                eph = np.append(eph, self.e) 

                self.T_GD_1 = float(self.rename_exp_designator(rinex_str[42:61]))
                eph = np.append(eph, self.T_GD_1)

                self.T_GD_2 = float(self.rename_exp_designator(rinex_str[61:80]))
                eph = np.append(eph, self.T_GD_2)

                # BROADCAST ORBIT-7

                rinex_str = file.read(max_char)
                num_rinex_str += 1

                if (rinex_str == None):
                    break

                self.SOW = float(self.rename_exp_designator(rinex_str[4:23]))
                eph = np.append(eph, self.SOW)

                self.AODC = float(self.rename_exp_designator(rinex_str[24:42]))
                eph = np.append(eph, self.AODC)

                health_flag = 1
                eph = np.append(eph, health_flag)

                NUM_SAT += 1

                if all_eph.size != 0:
                    mode = 0

            else:
                num_rinex_str += 7
            