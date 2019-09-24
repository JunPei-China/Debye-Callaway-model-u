# 质量波动和应力波动模型
import os,yaml,sys
import numpy as np
import scipy.constants as C


# 定义全局变量
pi = C.pi
h = C.Planck
hbar = h/(2*pi)
R = C.R
k = C.k

prompt = ">>>"

def periodic_table():
    CurrentPath=os.getcwd()
    YamlFile=os.path.join(CurrentPath,"periodic-table.yaml")

    with open(YamlFile,"r") as f:
        PeriodicTable = yaml.load(f,Loader=yaml.FullLoader)
    return PeriodicTable

def measure_data():
    LogitudinalSoundVelocity = float(input("请输入纵波声速(单位:m/s)"))
    TransverseSoundVelocity = float(input("请输入横波声速(单位:m/s)"))
    Density = float(input("请输入样品的密度(单位:g/cm^-3)"))
    NumberAtoms = float(input("请输入单胞中原子的个数"))
    VolumeCell = float(input("请输入单胞的体积(单位:A^-3)"))
    AverageAtomicVolume = VolumeCell/NumberAtoms
    CurrentPath = os.getcwd()
    AllFiles = os.listdir(CurrentPath)
    DataFile = []
    FileNum = 0
    for file in AllFiles:
        if os.path.splitext(file)[-1] == ".dat":
            FileNum += 1
            print(file)
    if FileNum == 0:
        print("请准备基体晶格热导率与温度的数据文件。格式:xxx.dat")
        try:
            sys.exit(0)
        except:
            print("die")
        finally:
            print("clean up.")
    print("请选择基体晶格热导率与温度的数据文件。格式:xxx.dat")
    MeasuredFile = os.path.join(CurrentPath,str(input(prompt)))
       
    LatticeThermalconductivity = np.loadtxt(MeasuredFile)


    # 数据规范化
    Measure = [LogitudinalSoundVelocity,TransverseSoundVelocity,Density,NumberAtoms,VolumeCell,AverageAtomicVolume,LatticeThermalconductivity]
    return Measure

def intrinsic_data():
    PeriodicTable = periodic_table()
    NumberAtom = int(input("请输入晶格中不同原子位置的种类个数")) 
    Elements = []
    for x in range(NumberAtom):
        MatrixElements = []
        NumberElement = int(input("请输入第%d种位置中元素种类个数"%(x+1)))
        for y in range(NumberElement):    
            NameElements = str(input("请输入第%d个位置第%d种元素的元素符号"%(x+1,y+1)))
            ElementSubscript = float(input("请输入第%d个位置第%d种元素的理论下标量"%(x+1,y+1)))
            AtomicMass = PeriodicTable[NameElements]["RelativeAtomicMass"]
            AtomicRadius = PeriodicTable[NameElements]["AtomicRadius"]["value"]
            MatrixElements.append([NameElements,ElementSubscript,AtomicMass,AtomicRadius])
        ElementMole = np.sum(np.array([i[1] for i in MatrixElements]))            
        ElementMass = np.dot(np.array([i[2] for i in MatrixElements]),np.array([i[1] for i in MatrixElements]))                

        AverageElementMass = ElementMass/ElementMole
        AverageMoleRadius = np.dot(AtomicRadius,ElementSubscript)/ElementMole        

        DopingOrNot = str(input("该位置是否有掺杂？y/n"))
        Doping = []
        if DopingOrNot[0].lower() == "y":
            NameDoping = str(input("请输入第%d个位置掺杂元素的元素符号"%(x+1)))
            DopingSubscript = float(input("请输入第%d个位置掺杂元素的下标量"%(x+1))) 
            DopingMass = PeriodicTable[NameDoping]["RelativeAtomicMass"]
            DopingRadius = PeriodicTable[NameDoping]["AtomicRadius"]["value"]
            Doping = [1,DopingSubscript,DopingMass,DopingRadius]
        else:
            Doping = [0,0,0,0]

        Matrix = [x,ElementMole,AverageElementMass,AverageMoleRadius]
        Matrix.extend(Doping)
    
        Elements.append(Matrix)
    print(Elements)

    return Elements

# Information = [x,ElementMole,AverageElementMass,AverageMoleRadius,\NameDoping,DopingSubscript,DopingMass,DopingRadius]


def average_sound_velocity(data):
    # data格式为: measure_data() Measure=[LogitudinalSoundVelocity,TransverseSoundVelocity,Density,\ NumberAtoms,VolumeCell,AverageAtomicVolume,LatticeThermalconductivity]
    LogitudinalSoundVelocity = data[0]
    TransverseSoundVelocity = data[1]
    
    AverageSoundVelocity = np.power(1/3*(1/np.power(LogitudinalSoundVelocity,3)+2/np.power(TransverseSoundVelocity,3)),-1/3)

    return AverageSoundVelocity


def Debye_temperature(data):
    # data格式为: measure_data() Measure=[LogitudinalSoundVelocity,
    # TransverseSoundVelocity,DensityNumberAtoms,VolumeCell,
    # AverageAtomicVolume,LatticeThermalconductivity]
    Density = data[2]
    NumberAtoms = data[3]
    VolumeCell = data[4]
    AverageSoundVelocity = average_sound_velocity(data)

    DebyeTemperature = h / k* np.power(3*NumberAtoms / (4*pi*VolumeCell),1/3) * AverageSoundVelocity
    
    return DebyeTemperature*1E10


def shear_modulus(data):
    # data格式为: measure_data() Measure=[LogitudinalSoundVelocity,\
    # TransverseSoundVelocity,Density,NumberAtoms,VolumeCell,\
    # AverageAtomicVolume,LatticeThermalconductivity]
    TransverseSoundVelocity = data[1]
    Density = data[2]
    
    ShearModulus = Density*np.power(TransverseSoundVelocity,2)
    
    return ShearModulus*1E-6

def bulk_modulus(data):
    # data格式为: measure_data() Measure=[LogitudinalSoundVelocity,
    #TransverseSoundVelocity,Density,NumberAtoms,VolumeCell,
    # AverageAtomicVolume,LatticeThermalconductivity]
    LogitudinalSoundVelocity = data[0]
    TransverseSoundVelocity = data[1]
    Density = data[2]
    
    BulkModulus = Density*(np.power(LogitudinalSoundVelocity,2)-4/3*np.power(TransverseSoundVelocity,2))
    return BulkModulus*1E-6

def Young_modulus(data):
    # data格式为: measure_data() Measure=[LogitudinalSoundVelocity,
    # TransverseSoundVelocity,Density,NumberAtoms,VolumeCell,
    # AverageAtomicVolume,LatticeThermalconductivity]
    LogitudinalSoundVelocity = data[0]
    TransverseSoundVelocity = data[1]
    Density = data[2]
    ShearModulus = shear_modulus(data)
    
    YoungModulus = ShearModulus*(3*np.power(LogitudinalSoundVelocity,2)-np.power(TransverseSoundVelocity,2))/(np.power(LogitudinalSoundVelocity,2)-np.power(TransverseSoundVelocity,2))
    
    return YoungModulus*1E-6

def Poisson_ratio(data):
    # data格式为: measure_data() Measure=[LogitudinalSoundVelocity,
    # TransverseSoundVelocity,Density,NumberAtoms,VolumeCell,
    # AverageAtomicVolume,LatticeThermalconductivity]
    LogitudinalSoundVelocity = data[0]
    TransverseSoundVelocity = data[1]
    Ratio = LogitudinalSoundVelocity/TransverseSoundVelocity
    
    PoissonRatio = (np.power(Ratio,2)-2)/(2*(np.power(Ratio,2)-1))
    
    return PoissonRatio

def K_modulus(data):
    # data格式为: measure_data() Measure=[LogitudinalSoundVelocity,
    # TransverseSoundVelocity,Density,NumberAtoms,VolumeCell,
    # AverageAtomicVolume,LatticeThermalconductivity]
    YoungModulus = Young_modulus(data)
    PoissonRatio = Poisson_ratio(data)
    
    KModulus = YoungModulus/(3*(1-2*PoissonRatio))
    
    return KModulus

def Gruneisen_parameter(data):
    # data格式为: measure_data() Measure=[LogitudinalSoundVelocity,
    # TransverseSoundVelocity,Density,NumberAtoms,VolumeCell,
    # AverageAtomicVolume,LatticeThermalconductivity]
    
    PoissonRatio = Poisson_ratio(data)

    
    GruneisenParameter = 3/2*((1+PoissonRatio)/(2-3*PoissonRatio))
    
    return GruneisenParameter

def strain_field_related_adjustable_parameter(data):
    # data格式为: measure_data() Measure=[LogitudinalSoundVelocity,
    # TransverseSoundVelocity,Density,NumberAtoms,VolumeCell,
    # AverageAtomicVolume,LatticeThermalconductivity]
    PoissonRatio = Poisson_ratio(data)
    GruneisenParameter = Gruneisen_parameter(data)
    Varepsilon = 2/9*np.power((6.4*GruneisenParameter*(1+PoissonRatio))/(1-PoissonRatio),2) 
    return Varepsilon

def thermal_conductivity_glass_limit_Cahill(data):
    # data格式为: measure_data() Measure=[LogitudinalSoundVelocity,
    # TransverseSoundVelocity,Density,NumberAtoms,VolumeCell,
    # AverageAtomicVolume,LatticeThermalconductivity]
    LogitudinalSoundVelocity = data[0]
    TransverseSoundVelocity = data[1]
    NumberAtoms = data[3]
    VolumeCell = data[4]
    
    CahillModel = 1/2*np.power(pi/6,1/3)*k*np.power(VolumeCell/NumberAtoms,-2/3)*(2*TransverseSoundVelocity+LogitudinalSoundVelocity)
    
    return CahillModel*1E20


def thermal_conductivity_glass_limit_Tan(data):
    # data格式为: measure_data() Measure=[LogitudinalSoundVelocity,
    # TransverseSoundVelocity,Density,NumberAtoms,VolumeCell,
    # AverageAtomicVolume,LatticeThermalconductivity]
    LogitudinalSoundVelocity = data[0]
    TransverseSoundVelocity = data[1]
    NumberAtoms = data[3]
    VolumeCell = data[4]
    
    AverageSoundVelocity = average_sound_velocity(data)
    TanModel = pi/4*k*np.power(VolumeCell/NumberAtoms,-2/3)*AverageSoundVelocity
    
    return TanModel*1E20

def Gamma_mass_and_strain(data1,data2):
    # data格式为: measure_data() Measure=[LogitudinalSoundVelocity,
    # TransverseSoundVelocity,Density,NumberAtoms,VolumeCell,
    # AverageAtomicVolume,LatticeThermalconductivity]

    # data2 格式为: intrinsic_data() Information = [x,ElementMole,
    # AverageElementMass,AverageMoleRadius,NameDoping,DopingSubscript,
    # DopingMass,DopingRadius]
    Varepsilon = strain_field_related_adjustable_parameter(data1)
    Information = np.array(data2)

    TotalMole = np.sum(Information[:,1])

    AverageMassAtom = ((Information[:,1]-Information[:,5])*Information[:,2]+Information[:,5]*Information[:,6])/Information[:,1]

    AverageMassElement = np.sum(AverageMassAtom*Information[:,1])/TotalMole

    AverageradiusAtom = ((Information[:,1]-Information[:,5])*Information[:,3]+Information[:,5]*Information[:,7])/Information[:,1]    
                        
    GammaMass = np.sum(Information[:,1]*np.power(AverageMassAtom/AverageMassElement,2)*(Information[:,1]-Information[:,5])/Information[:,1]*Information[:,5]/Information[:,1]*np.power((Information[:,2]-Information[:,6])/AverageMassAtom,2))/TotalMole
    GammaStrain = np.sum(Information[:,1]*np.power(AverageMassAtom/AverageMassElement,2)*(Information[:,1]-Information[:,5])/Information[:,1]*Information[:,5]/Information[:,1]*Varepsilon*np.power((Information[:,3]-Information[:,7])/AverageradiusAtom,2))/TotalMole

    GammaTotal = GammaMass+GammaStrain
    return (GammaMass,GammaStrain,GammaTotal)

def disorder_scattering_parameter(data1,data2):
    # data格式为: measure_data() Measure=[LogitudinalSoundVelocity,
    # TransverseSoundVelocity,Density,NumberAtoms,VolumeCell,
    # AverageAtomicVolume,LatticeThermalconductivity]

    # data2 格式为: intrinsic_data() Information = [x,ElementMole,
    # AverageElementMass,AverageMoleRadius,NameDoping,DopingSubscript,
    # DopingMass,DopingRadius]

    Measure = np.array(data1)
    Information = np.array(data2)
    DebyeTemperature = Debye_temperature(Measure)
    AverageAtomicVolume = Measure[5]*1E-30
    AverageSoundVelocity = average_sound_velocity(Measure)
    GammaTotal = Gamma_mass_and_strain(Measure,Information)[2]
    LatticeThermalconductivity = Measure[6][:,1]
    Temperature = Measure[6][:,0]
    u = np.sqrt((pi*pi*DebyeTemperature*AverageAtomicVolume)/(h*AverageSoundVelocity*AverageSoundVelocity)*LatticeThermalconductivity*GammaTotal)
    
    return np.array([Temperature,u]).T

def lattice_thermal_conductivity_callaway_model(data1,data2):
    # data格式为: measure_data() Measure=[LogitudinalSoundVelocity,
    # TransverseSoundVelocity,Density,NumberAtoms,VolumeCell,
    # AverageAtomicVolume,LatticeThermalconductivity]

    # data2 格式为: intrinsic_data() Information = [x,ElementMole,
    # AverageElementMass,AverageMoleRadius,NameDoping,DopingSubscript,
    # DopingMass,DopingRadius]
    Measure = np.array(data1)
    Information = np.array(data2)
    PureLatticeThermalconductivity = Measure[6][:,1]
    Temperature = Measure[6][:,0]
    u = disorder_scattering_parameter(Measure,Information)[:,1]

    LatticeThermalconductivity = np.arctan(u)*PureLatticeThermalconductivity/u
    
    return np.array([Temperature,LatticeThermalconductivity]).T

def test():
    print("*"*40+"参数输入"+"*"*40)
    print("")    
    Measure = measure_data()
    Information = intrinsic_data()
    
    print("")
    print("")
    print("")
    print("*"*40+"结果输出"+"*"*40)
    AverageSoundVelocity = average_sound_velocity(Measure)    
    DebyeTemperature = Debye_temperature(Measure)
    ShearModulus = shear_modulus(Measure)
    YoungModulus = Young_modulus(Measure)
    BulkModulus = bulk_modulus(Measure)
    KModulus = K_modulus(Measure)
    PoissonRatio = Poisson_ratio(Measure)
    GruneisenParameter = Gruneisen_parameter(Measure)
    Varepsilon = strain_field_related_adjustable_parameter(Measure)
    ThermalConductivityGlassLimitCahill = thermal_conductivity_glass_limit_Cahill(Measure)
    ThermalConductivityGlassLimitTan = thermal_conductivity_glass_limit_Tan(Measure)
    
    #输出声速测试结果
    print("")
    print("平均声速为: %f m/s"%AverageSoundVelocity)
    print("德拜温度为: %f K"%DebyeTemperature)
    print("剪切模量G为: %f MPa"%ShearModulus)
    print("杨氏模量E为: %f MPa"%YoungModulus)
    print("块体模量B为: %f MPa"%BulkModulus)
    print("体积模量K为: %f MPa"%KModulus)
    print("泊松比: %f"%PoissonRatio)
    print("格林艾森常数: %f"%GruneisenParameter)
    print("应力波动调整参数varepsilon: %f"%Varepsilon)
    print("晶格热导非晶极限(Cahill模型): %f W/m/K"%ThermalConductivityGlassLimitCahill)
    print("晶格热导非晶极限(Tan文章模型): %f W/m/K"%ThermalConductivityGlassLimitTan)
    
    #输出应力波动和质量波动
    GammaMass,GammaStrain,GammaTotal = Gamma_mass_and_strain(Measure,Information)
    print("")
    print("质量波动:%f"%GammaMass)
    print("应力波动： %f"%GammaStrain)
    print("总波动: %f"%GammaTotal)
    
    #输出无序度因子
    u = disorder_scattering_parameter(Measure,Information)    
    print("")
    print("无序度因子")
    print("T-K u")
    print(u)

    
    LatticeThermalconductivity = lattice_thermal_conductivity_callaway_model(Measure,Information)    
    
    # 输出
    print("")
    print("固溶体的晶格热导率(Callaway模型)为:")
    print("T-K k-W/m/K")
    print(LatticeThermalconductivity)
    ExitInput=input("按任意键退出")
                            
if __name__ == "__main__":
    test()
    

