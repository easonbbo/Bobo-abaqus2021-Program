C. Fortran 语言中实现导向辊运动控制的子程序 
C          user amplitude subroutine 
      
      !子程序标准开头与参数定义
      
      Subroutine VUAMP( 
C                    c 
     *          ampName,time,ampValueOld,dt, nprops,props,nSvars,svars,!svars，用于在不同的增量步之间保存数据，让程序记住上一步发生了什么，实现连续控制。
     *          lFlagsInfo, 
     *          nSensor, sensorValues, sensorNames,   
     *          jSensorLookUpTable,   
c                    to be defined 
     *          ampValueNew,   
     *          lFlagsDefine, 
     *          AmpDerivative, AmpSecDerivative, AmpIncIntegral)             
      include 'vaba_param.inc' 
     !!幅值名称ampName（区分不同轧辊的控制通道）、时间time、
     !! 传感器数据sensorValues（对应论文的虚拟传感器）、
     !! 状态变量svars（保存轧制历史数据）
      
      !输出参数ampValueNew是计算出的轧辊运动速度，ABAQUS 会用该值驱动轧辊运动，实现自适应控制。
      
      
      
      !数组与变量声明
C          svars - additional state variables, similar to (V)UEL 
      dimension sensorValues(nSensor), svars(nSvars) !前一个长度 = 传感器总数 nSensor,后一个长度 = 状态变量总数 nSvars
      character*80 sensorNames(nSensor) !长度80个字符
      character*80 ampName !你要控制的设备名,A1 = 导向辊 Y 方向
	
C          time indices  
      parameter( iStepTime   = 1, 
     *           iTotalTime  = 2, 
     *           nTime       = 2) 
C          flags passed in for information 
      parameter( iInitialization = 1, 
     *           iRegularInc     = 2, 
     *           nFlagsInfo      = 2) 
C          optional flags to be defined 
      parameter(  iComputeDeriv    = 1, 
     *            iComputeSecDeriv = 2, 
     *            iComputeInteg    = 3, 
     *            iStopAnalysis    = 4, 
     *            iConcludeStep    = 5, 
     *            nFlagsDefine     = 5) 
      
      dimension time(nTime), lFlagsInfo(nFlagsInfo), 
     *                    lFlagsDefine(nFlagsDefine) 
 
      dimension jSensorLookUpTable(*) 
c     dimension x(63),y(63) 
	real(8)::x(200),y(200)  
      !双精度数组，用于存储环件外表面 8 个截面、共 192 个节点的 X/Y 坐标                                                                                 
      real(8)         t,
     *                curve,
     *                Drivex,
     *                Corex,
     *                Cupy,
     *                Cdowny,
     *                Cupx,
     *                Cdownx,   
     *                Rideal,
     *                Vd,        !环件长大速度
     *                Vd0,        
     *                Vf,        !芯辊径向
     *                Va,        !锥辊轴向
     *                bdet,      !瞬时壁厚测量值
     *                b,         !瞬时壁厚         
     *                h,         !瞬时高度
     *                L, 
     *                R,                                                                                                      
     *              angle,   
     *                w ,
     *            Lcone,
     *			GRotRigD,
     *                 n,
     *                 lt,   
     *                 Rn,
     *                 Rt
      
	
      parameter ( 
     *            Acone=16/180*3.141592654,   !锥辊锥角,用于锥辊运动的速度分解
     *            Vdrive=700,                 !驱动辊转速
     *            Ridle=150,                  !芯辊半径                                                                        
     *            Rguide=100,                 !导辊半径                                                 
     *            Rmain=500,                  !驱动辊半径                                                                                                                                                                                                                              
     *			GRotXo=-600,                !与左导向辊相关的铰链位置
     *			GRotYo=-600,
      
     *            Guide_X=678,              !L导辊X坐标
     *            Guide_Z=-362.13,            !L导辊Z坐标,
        !导向辊的旋转铰链中心坐标，Rg为导向臂长度
      
     *	        Rg=sqrt((GRotXo-Guide_X)**2+(GRotYo-Guide_Z)**2),    
      
     *            ti=0.05, 
      !子程序每 0.05s 和 ABAQUS 主程序交互一次，采集传感器数据并更新轧辊速度
     *			pi=3.1415926,
     *           ! Vblank=3.9740594,
     *			b0=56,     !初始壁厚
     *			bf=25,       !结束壁厚
     *			h0=930,      !初始高度
     *			hf=930,      !初始高度
     *			R0=277,      !初始外半径,
     *			Rf=470.5,    !结束半径,大端
     *            bi=218.75,
     *            hi=226.)
     
      
      !仿真初始化步处理
      lFlagsDefine(iComputeDeriv)    = 1 
      lFlagsDefine(iComputeSecDeriv) = 1  
      !开启幅值的一阶、二阶导数计算，保证轧辊运动的平滑性，避免冲击
      
      if(lFlagsInfo(iInitialization).eq.1) then 
c NOTE: the model defines the amplitudes with variables=20, so svars
c has only 20 slots. Per-channel clocks must stay inside 1..20:
c A1->1, A2->2, A3->3, A4->4, A3-1->12, A4-1->13, A5->14.
      !当仿真进入初始化步时，对状态变量svars赋初始值
			select case(ampName)
      !区分不同轧辊的控制通道，给导向辊、锥辊赋初始坐标 / 转速，对应轧制前的设备初始位置设置
				case('A1') 
				svars(5)=-362.13     !导辊Y坐标(左边的导辊)
				svars(1)=-ti
				case('A2') 
				svars(6)=678         !导辊X坐标
				svars(2)=-ti
                  case('A3')
				svars(7)=2*R0
				svars(3)=-ti
                  case('A4') 
				svars(8)=2*R0
				svars(4)=-ti
                  case('A3-1') 
				svars(10)=2*R0
				svars(12)=-ti
                  case('A4-1') 
				svars(11)=2*R0
				svars(13)=-ti
                  case('A5') 
                  ampValueNew=6.0299
                  svars(14)=-ti
			end select           
	else 
c Each amplitude gets its own control clock. VUAMP is called once per
c amplitude per increment; sharing one clock means only the first
c channel updates and all others keep their init value (0, or 6.0299
c for A5): cone rollers then spin but never translate.
c Clock slots (must stay within nSvars=20, see *Amplitude variables):
c A1->1, A2->2, A3->3, A4->4, A3-1->12, A4-1->13, A5->14.
				kTick=0
				if(ampName(1:2).eq.'A1') then
					kTick=1
				else if(ampName(1:2).eq.'A2') then
					kTick=2
				else if(ampName(1:3).eq.'A3-') then
					kTick=12
				else if(ampName(1:2).eq.'A3') then
					kTick=3
				else if(ampName(1:3).eq.'A4-') then
					kTick=13
				else if(ampName(1:2).eq.'A4') then
					kTick=4
				else
					kTick=14
				end if
				if(kTick.gt.0 .and. (svars(kTick).le.0.d0 .or.
     *				time(iTotalTime)-svars(kTick).ge.ti-1.d-6)) then
c The initialization phase of VUAMP is only called once, with the first
c amplitude name, so every channel must self-initialize on its own first
c regular call. svars(kTick)<=0 marks "first call" for this channel.
						lFirst=0
						if(svars(kTick).le.0.d0) lFirst=1
						t=time(iTotalTime)-svars(kTick)
						if(t.le.0.d0) t=ti
						svars(kTick)=time(iTotalTime)
      !控制周期判断与时间更新
                          
                          
c          get sensor value 计算R
	x(1)    = vGetSensorValue('XA1', 
     * jSensorLookUpTable, 
     * sensorValues)
      !jSensorLookUpTable：整数数组，传感器名称到数值索引的快速查找表，
      !配合 vGetSensorValue 函数使用，通过传感器名字直接找到它在 sensorValues 中的值
	x(2)    = vGetSensorValue('XA2', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(3)    = vGetSensorValue('XA3', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(4)    = vGetSensorValue('XA4', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(5)    = vGetSensorValue('XA5', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(6)= vGetSensorValue('XA6', 
     * jSensorLookUpTable,  
     * sensorValues) 
	x(7)    = vGetSensorValue('XA7', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(8)    = vGetSensorValue('XA8', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(9)    = vGetSensorValue('XA9', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(10)    = vGetSensorValue('XA10', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(11)    = vGetSensorValue('XA11', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(12)= vGetSensorValue('XA12', 
     * jSensorLookUpTable,  
     * sensorValues)
	x(13)    = vGetSensorValue('XA13', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(14)    = vGetSensorValue('XA14', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(15)    = vGetSensorValue('XA15', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(16)    = vGetSensorValue('XA16', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(17)    = vGetSensorValue('XA17', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(18)= vGetSensorValue('XA18', 
     * jSensorLookUpTable,  
     * sensorValues) 
	x(19)    = vGetSensorValue('XA19', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(20)    = vGetSensorValue('XA20', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(21)    = vGetSensorValue('XA21', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(22)    = vGetSensorValue('XA22', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(23)    = vGetSensorValue('XA23', 
     *jSensorLookUpTable, 
     *sensorValues) 
     	x(24)= vGetSensorValue('XA24', 
     *jSensorLookUpTable, 
     * sensorValues)
      
      x(25)    = vGetSensorValue('XB1', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(26)    = vGetSensorValue('XB2', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(27)    = vGetSensorValue('XB3', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(28)    = vGetSensorValue('XB4', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(29)    = vGetSensorValue('XB5', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(30)= vGetSensorValue('XB6', 
     * jSensorLookUpTable,  
     * sensorValues) 
	x(31)    = vGetSensorValue('XB7', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(32)    = vGetSensorValue('XB8', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(33)    = vGetSensorValue('XB9', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(34)    = vGetSensorValue('XB10', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(35)    = vGetSensorValue('XB11', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(36)= vGetSensorValue('XB12', 
     * jSensorLookUpTable,  
     * sensorValues)
	x(37)    = vGetSensorValue('XB13', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(38)    = vGetSensorValue('XB14', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(39)    = vGetSensorValue('XB15', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(40)    = vGetSensorValue('XB16', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(41)    = vGetSensorValue('XB17', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(42)= vGetSensorValue('XB18', 
     * jSensorLookUpTable,  
     * sensorValues) 
	x(43)    = vGetSensorValue('XB19', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(44)    = vGetSensorValue('XB20', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(45)    = vGetSensorValue('XB21', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(46)    = vGetSensorValue('XB22', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(47)    = vGetSensorValue('XB23', 
     *jSensorLookUpTable, 
     *sensorValues) 
     	x(48)= vGetSensorValue('XB24', 
     *jSensorLookUpTable, 
     * sensorValues)
      
       x(49)  = vGetSensorValue('XC1', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(50)    = vGetSensorValue('XC2', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(51)    = vGetSensorValue('XC3', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(52)    = vGetSensorValue('XC4', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(53)    = vGetSensorValue('XC5', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(54)= vGetSensorValue('XC6', 
     * jSensorLookUpTable,  
     * sensorValues) 
	x(55)    = vGetSensorValue('XC7', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(56)    = vGetSensorValue('XC8', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(57)    = vGetSensorValue('XC9', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(58)    = vGetSensorValue('XC10', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(59)    = vGetSensorValue('XC11', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(60)= vGetSensorValue('XC12', 
     * jSensorLookUpTable,  
     * sensorValues)
	x(61)    = vGetSensorValue('XC13', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(62)    = vGetSensorValue('XC14', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(63)    = vGetSensorValue('XC15', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(64)    = vGetSensorValue('XC16', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(65)    = vGetSensorValue('XC17', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(66)= vGetSensorValue('XC18', 
     * jSensorLookUpTable,  
     * sensorValues) 
	x(67)    = vGetSensorValue('XC19', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(68)    = vGetSensorValue('XC20', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(69)    = vGetSensorValue('XC21', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(70)    = vGetSensorValue('XC22', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(71)    = vGetSensorValue('XC23', 
     *jSensorLookUpTable, 
     *sensorValues) 
     	x(72)= vGetSensorValue('XC24', 
     *jSensorLookUpTable, 
     * sensorValues)
      
      x(73)    = vGetSensorValue('XD1', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(74)    = vGetSensorValue('XD2', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(75)    = vGetSensorValue('XD3', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(76)    = vGetSensorValue('XD4', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(77)    = vGetSensorValue('XD5', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(78)= vGetSensorValue('XD6', 
     * jSensorLookUpTable,  
     * sensorValues) 
	x(79)    = vGetSensorValue('XD7', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(80)    = vGetSensorValue('XD8', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(81)    = vGetSensorValue('XD9', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(82)    = vGetSensorValue('XD10', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(83)    = vGetSensorValue('XD11', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(84)= vGetSensorValue('XD12', 
     * jSensorLookUpTable,  
     * sensorValues)
	x(85)    = vGetSensorValue('XD13', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(86)    = vGetSensorValue('XD14', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(87)    = vGetSensorValue('XD15', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(88)    = vGetSensorValue('XD16', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(89)    = vGetSensorValue('XD17', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(90)= vGetSensorValue('XD18', 
     * jSensorLookUpTable,  
     * sensorValues) 
	x(91)    = vGetSensorValue('XD19', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(92)    = vGetSensorValue('XD20', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(93)    = vGetSensorValue('XD21', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(94)    = vGetSensorValue('XD22', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(95)    = vGetSensorValue('XD23', 
     *jSensorLookUpTable, 
     *sensorValues) 
     	x(96)= vGetSensorValue('XD24', 
     *jSensorLookUpTable, 
     * sensorValues)
      
      x(97)    = vGetSensorValue('XE1', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(98)    = vGetSensorValue('XE2', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(99)    = vGetSensorValue('XE3', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(100)    = vGetSensorValue('XE4', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(101)    = vGetSensorValue('XE5', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(102)= vGetSensorValue('XE6', 
     * jSensorLookUpTable,  
     * sensorValues) 
	x(103)    = vGetSensorValue('XE7', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(104)    = vGetSensorValue('XE8', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(105)    = vGetSensorValue('XE9', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(106)    = vGetSensorValue('XE10', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(107)    = vGetSensorValue('XE11', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(108)= vGetSensorValue('XE12', 
     * jSensorLookUpTable,  
     * sensorValues)
	x(109)    = vGetSensorValue('XE13', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(110)    = vGetSensorValue('XE14', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(111)    = vGetSensorValue('XE15', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(112)    = vGetSensorValue('XE16', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(113)    = vGetSensorValue('XE17', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(114)= vGetSensorValue('XE18', 
     * jSensorLookUpTable,  
     * sensorValues) 
	x(115)    = vGetSensorValue('XE19', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(116)    = vGetSensorValue('XE20', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(117)    = vGetSensorValue('XE21', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(118)    = vGetSensorValue('XE22', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(119)    = vGetSensorValue('XE23', 
     *jSensorLookUpTable, 
     *sensorValues) 
     	x(120)= vGetSensorValue('XE24', 
     *jSensorLookUpTable, 
     * sensorValues)
      
      x(121)    = vGetSensorValue('XF1', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(122)    = vGetSensorValue('XF2', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(123)    = vGetSensorValue('XF3', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(124)    = vGetSensorValue('XF4', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(125)    = vGetSensorValue('XF5', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(126)= vGetSensorValue('XF6', 
     * jSensorLookUpTable,  
     * sensorValues) 
	x(127)    = vGetSensorValue('XF7', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(128)    = vGetSensorValue('XF8', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(129)    = vGetSensorValue('XF9', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(130)    = vGetSensorValue('XF10', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(131)    = vGetSensorValue('XF11', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(132)= vGetSensorValue('XF12', 
     * jSensorLookUpTable,  
     * sensorValues)
	x(133)    = vGetSensorValue('XF13', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(134)    = vGetSensorValue('XF14', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(135)    = vGetSensorValue('XF15', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(136)    = vGetSensorValue('XF16', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(137)    = vGetSensorValue('XF17', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(138)= vGetSensorValue('XF18', 
     * jSensorLookUpTable,  
     * sensorValues) 
	x(139)    = vGetSensorValue('XF19', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(140)    = vGetSensorValue('XF20', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(141)    = vGetSensorValue('XF21', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(142)    = vGetSensorValue('XF22', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(143)    = vGetSensorValue('XF23', 
     *jSensorLookUpTable, 
     *sensorValues) 
     	x(144)= vGetSensorValue('XF24', 
     *jSensorLookUpTable, 
     * sensorValues)
      
       x(145)    = vGetSensorValue('XG1', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(146)    = vGetSensorValue('XG2', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(147)    = vGetSensorValue('XG3', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(148)    = vGetSensorValue('XG4', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(149)    = vGetSensorValue('XG5', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(150)= vGetSensorValue('XG6', 
     * jSensorLookUpTable,  
     * sensorValues) 
	x(151)    = vGetSensorValue('XG7', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(152)    = vGetSensorValue('XG8', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(153)    = vGetSensorValue('XG9', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(154)    = vGetSensorValue('XG10', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(155)    = vGetSensorValue('XG11', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(156)= vGetSensorValue('XG12', 
     * jSensorLookUpTable,  
     * sensorValues)
	x(157)    = vGetSensorValue('XG13', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(158)    = vGetSensorValue('XG14', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(159)    = vGetSensorValue('XG15', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(160)    = vGetSensorValue('XG16', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(161)    = vGetSensorValue('XG17', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(162)= vGetSensorValue('XG18', 
     * jSensorLookUpTable,  
     * sensorValues) 
	x(163)    = vGetSensorValue('XG19', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(164)    = vGetSensorValue('XG20', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(165)    = vGetSensorValue('XG21', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(166)    = vGetSensorValue('XG22', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(167)    = vGetSensorValue('XG23', 
     *jSensorLookUpTable, 
     *sensorValues) 
     	x(168)= vGetSensorValue('XG24', 
     *jSensorLookUpTable, 
     * sensorValues)
      
      x(169)    = vGetSensorValue('XH1', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(170)    = vGetSensorValue('XH2', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(171)    = vGetSensorValue('XH3', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(172)    = vGetSensorValue('XH4', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(173)    = vGetSensorValue('XH5', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(174)= vGetSensorValue('XH6', 
     * jSensorLookUpTable,  
     * sensorValues) 
	x(175)    = vGetSensorValue('XH7', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(176)    = vGetSensorValue('XH8', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(177)    = vGetSensorValue('XH9', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(178)    = vGetSensorValue('XH10', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(179)    = vGetSensorValue('XH11', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(180)= vGetSensorValue('XH12', 
     * jSensorLookUpTable,  
     * sensorValues)
	x(181)    = vGetSensorValue('XH13', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(182)    = vGetSensorValue('XH14', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(183)    = vGetSensorValue('XH15', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(184)    = vGetSensorValue('XH16', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(185)    = vGetSensorValue('XH17', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(186)= vGetSensorValue('XH18', 
     * jSensorLookUpTable,  
     * sensorValues) 
	x(187)    = vGetSensorValue('XH19', 
     * jSensorLookUpTable, 
     * sensorValues) 
	x(188)    = vGetSensorValue('XH20', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	x(189)    = vGetSensorValue('XH21', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(190)    = vGetSensorValue('XH22', 
     *jSensorLookUpTable, 
     *sensorValues) 
	x(191)    = vGetSensorValue('XH23', 
     *jSensorLookUpTable, 
     *sensorValues) 
     	x(192)= vGetSensorValue('XH24', 
     *jSensorLookUpTable, 
     * sensorValues)
      
	
      
	y(1)    = vGetSensorValue('YA1', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(2)    = vGetSensorValue('YA2', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(3)    = vGetSensorValue('YA3', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(4)    = vGetSensorValue('YA4', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(5)    = vGetSensorValue('YA5', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(6)= vGetSensorValue('YA6', 
     * jSensorLookUpTable,  
     * sensorValues) 
	y(7)    = vGetSensorValue('YA7', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(8)    = vGetSensorValue('YA8', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(9)    = vGetSensorValue('YA9', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(10)    = vGetSensorValue('YA10', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(11)    = vGetSensorValue('YA11', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(12)= vGetSensorValue('YA12', 
     * jSensorLookUpTable,  
     * sensorValues)
	y(13)    = vGetSensorValue('YA13', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(14)    = vGetSensorValue('YA14', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(15)    = vGetSensorValue('YA15', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(16)    = vGetSensorValue('YA16', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(17)    = vGetSensorValue('YA17', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(18)= vGetSensorValue('YA18', 
     * jSensorLookUpTable,  
     * sensorValues) 
	y(19)    = vGetSensorValue('YA19', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(20)    = vGetSensorValue('YA20', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(21)    = vGetSensorValue('YA21', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(22)    = vGetSensorValue('YA22', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(23)    = vGetSensorValue('YA23', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(24)= vGetSensorValue('YA24', 
     * jSensorLookUpTable,  
     * sensorValues)
     
      y(25)    = vGetSensorValue('YB1', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(26)    = vGetSensorValue('YB2', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(27)    = vGetSensorValue('YB3', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(28)    = vGetSensorValue('YB4', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(29)    = vGetSensorValue('YB5', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(30)= vGetSensorValue('YB6', 
     * jSensorLookUpTable,  
     * sensorValues) 
	y(31)    = vGetSensorValue('YB7', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(32)    = vGetSensorValue('YB8', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(33)    = vGetSensorValue('YB9', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(34)    = vGetSensorValue('YB10', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(35)    = vGetSensorValue('YB11', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(36)= vGetSensorValue('YB12', 
     * jSensorLookUpTable,  
     * sensorValues)
	y(37)    = vGetSensorValue('YB13', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(38)    = vGetSensorValue('YB14', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(39)    = vGetSensorValue('YB15', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(40)    = vGetSensorValue('YB16', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(41)    = vGetSensorValue('YB17', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(42)= vGetSensorValue('YB18', 
     * jSensorLookUpTable,  
     * sensorValues) 
	y(43)    = vGetSensorValue('YB19', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(44)    = vGetSensorValue('YB20', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(45)    = vGetSensorValue('YB21', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(46)    = vGetSensorValue('YB22', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(47)    = vGetSensorValue('YB23', 
     *jSensorLookUpTable, 
     *sensorValues) 
     	y(48)= vGetSensorValue('YB24', 
     *jSensorLookUpTable, 
     * sensorValues)   
 
      y(49)    = vGetSensorValue('YC1', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(50)    = vGetSensorValue('YC2', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(51)    = vGetSensorValue('YC3', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(52)    = vGetSensorValue('YC4', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(53)    = vGetSensorValue('YC5', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(54)= vGetSensorValue('YC6', 
     * jSensorLookUpTable,  
     * sensorValues) 
	y(55)    = vGetSensorValue('YC7', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(56)    = vGetSensorValue('YC8', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(57)    = vGetSensorValue('YC9', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(58)    = vGetSensorValue('YC10', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(59)    = vGetSensorValue('YC11', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(60)= vGetSensorValue('YC12', 
     * jSensorLookUpTable,  
     * sensorValues)
	y(61)    = vGetSensorValue('YC13', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(62)    = vGetSensorValue('YC14', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(63)    = vGetSensorValue('YC15', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(64)    = vGetSensorValue('YC16', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(65)    = vGetSensorValue('YC17', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(66)= vGetSensorValue('YC18', 
     * jSensorLookUpTable,  
     * sensorValues) 
	y(67)    = vGetSensorValue('YC19', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(68)    = vGetSensorValue('YC20', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(69)    = vGetSensorValue('YC21', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(70)    = vGetSensorValue('YC22', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(71)    = vGetSensorValue('YC23', 
     *jSensorLookUpTable, 
     *sensorValues) 
     	y(72)= vGetSensorValue('YC24', 
     *jSensorLookUpTable, 
     * sensorValues)
      
      y(73)    = vGetSensorValue('YD1', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(74)    = vGetSensorValue('YD2', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(75)    = vGetSensorValue('YD3', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(76)    = vGetSensorValue('YD4', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(77)    = vGetSensorValue('YD5', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(78)= vGetSensorValue('YD6', 
     * jSensorLookUpTable,  
     * sensorValues) 
	y(79)    = vGetSensorValue('YD7', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(80)    = vGetSensorValue('YD8', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(81)    = vGetSensorValue('YD9', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(82)    = vGetSensorValue('YD10', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(83)    = vGetSensorValue('YD11', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(84)= vGetSensorValue('YD12', 
     * jSensorLookUpTable,  
     * sensorValues)
	y(85)    = vGetSensorValue('YD13', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(86)    = vGetSensorValue('YD14', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(87)    = vGetSensorValue('YD15', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(88)    = vGetSensorValue('YD16', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(89)    = vGetSensorValue('YD17', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(90)= vGetSensorValue('YD18', 
     * jSensorLookUpTable,  
     * sensorValues) 
	y(91)    = vGetSensorValue('YD19', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(92)    = vGetSensorValue('YD20', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(93)    = vGetSensorValue('YD21', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(94)    = vGetSensorValue('YD22', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(95)    = vGetSensorValue('YD23', 
     *jSensorLookUpTable, 
     *sensorValues) 
     	y(96)= vGetSensorValue('YD24', 
     *jSensorLookUpTable, 
     * sensorValues)
      
      y(97)    = vGetSensorValue('YE1', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(98)    = vGetSensorValue('YE2', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(99)    = vGetSensorValue('YE3', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(100)    = vGetSensorValue('YE4', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(101)    = vGetSensorValue('YE5', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(102)= vGetSensorValue('YE6', 
     * jSensorLookUpTable,  
     * sensorValues) 
	y(103)    = vGetSensorValue('YE7', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(104)    = vGetSensorValue('YE8', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(105)    = vGetSensorValue('YE9', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(106)    = vGetSensorValue('YE10', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(107)    = vGetSensorValue('YE11', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(108)= vGetSensorValue('YE12', 
     * jSensorLookUpTable,  
     * sensorValues)
	y(109)    = vGetSensorValue('YE13', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(110)    = vGetSensorValue('YE14', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(111)    = vGetSensorValue('YE15', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(112)    = vGetSensorValue('YE16', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(113)    = vGetSensorValue('YE17', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(114)= vGetSensorValue('YE18', 
     * jSensorLookUpTable,  
     * sensorValues) 
	y(115)    = vGetSensorValue('YE19', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(116)    = vGetSensorValue('YE20', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(117)    = vGetSensorValue('YE21', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(118)    = vGetSensorValue('YE22', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(119)    = vGetSensorValue('YE23', 
     *jSensorLookUpTable, 
     *sensorValues) 
     	y(120)= vGetSensorValue('YE24', 
     *jSensorLookUpTable, 
     * sensorValues)
      
      y(121)    = vGetSensorValue('YF1', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(122)    = vGetSensorValue('YF2', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(123)    = vGetSensorValue('YF3', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(124)    = vGetSensorValue('YF4', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(125)    = vGetSensorValue('YF5', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(126)= vGetSensorValue('YF6', 
     * jSensorLookUpTable,  
     * sensorValues) 
	y(127)    = vGetSensorValue('YF7', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(128)    = vGetSensorValue('YF8', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(129)    = vGetSensorValue('YF9', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(130)    = vGetSensorValue('YF10', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(131)    = vGetSensorValue('YF11', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(132)= vGetSensorValue('YF12', 
     * jSensorLookUpTable,  
     * sensorValues)
	y(133)    = vGetSensorValue('YF13', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(134)    = vGetSensorValue('YF14', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(135)    = vGetSensorValue('YF15', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(136)    = vGetSensorValue('YF16', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(137)    = vGetSensorValue('YF17', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(138)= vGetSensorValue('YF18', 
     * jSensorLookUpTable,  
     * sensorValues) 
	y(139)    = vGetSensorValue('YF19', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(140)    = vGetSensorValue('YF20', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(141)    = vGetSensorValue('YF21', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(142)    = vGetSensorValue('YF22', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(143)    = vGetSensorValue('YF23', 
     *jSensorLookUpTable, 
     *sensorValues) 
     	y(144)= vGetSensorValue('YF24', 
     *jSensorLookUpTable, 
     * sensorValues)
      
       y(145)    = vGetSensorValue('YG1', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(146)    = vGetSensorValue('YG2', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(147)    = vGetSensorValue('YG3', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(148)    = vGetSensorValue('YG4', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(149)    = vGetSensorValue('YG5', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(150)= vGetSensorValue('YG6', 
     * jSensorLookUpTable,  
     * sensorValues) 
	y(151)    = vGetSensorValue('YG7', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(152)    = vGetSensorValue('YG8', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(153)    = vGetSensorValue('YG9', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(154)    = vGetSensorValue('YG10', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(155)    = vGetSensorValue('YG11', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(156)= vGetSensorValue('YG12', 
     * jSensorLookUpTable,  
     * sensorValues)
	y(157)    = vGetSensorValue('YG13', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(158)    = vGetSensorValue('YG14', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(159)    = vGetSensorValue('YG15', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(160)    = vGetSensorValue('YG16', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(161)    = vGetSensorValue('YG17', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(162)= vGetSensorValue('YG18', 
     * jSensorLookUpTable,  
     * sensorValues) 
	y(163)    = vGetSensorValue('YG19', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(164)    = vGetSensorValue('YG20', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(165)    = vGetSensorValue('YG21', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(166)    = vGetSensorValue('YG22', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(167)    = vGetSensorValue('YG23', 
     *jSensorLookUpTable, 
     *sensorValues) 
     	y(168)= vGetSensorValue('YG24', 
     *jSensorLookUpTable, 
     * sensorValues)
      
      y(169)    = vGetSensorValue('YH1', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(170)    = vGetSensorValue('YH2', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(171)    = vGetSensorValue('YH3', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(172)    = vGetSensorValue('YH4', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(173)    = vGetSensorValue('YH5', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(174)= vGetSensorValue('YH6', 
     * jSensorLookUpTable,  
     * sensorValues) 
	y(175)    = vGetSensorValue('YH7', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(176)    = vGetSensorValue('YH8', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(177)    = vGetSensorValue('YH9', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(178)    = vGetSensorValue('YH10', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(179)    = vGetSensorValue('YH11', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(180)= vGetSensorValue('YH12', 
     * jSensorLookUpTable,  
     * sensorValues)
	y(181)    = vGetSensorValue('YH13', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(182)    = vGetSensorValue('YH14', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(183)    = vGetSensorValue('YH15', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(184)    = vGetSensorValue('YH16', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(185)    = vGetSensorValue('YH17', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(186)= vGetSensorValue('YH18', 
     * jSensorLookUpTable,  
     * sensorValues) 
	y(187)    = vGetSensorValue('YH19', 
     * jSensorLookUpTable, 
     * sensorValues) 
	y(188)    = vGetSensorValue('YH20', 
     * jSensorLookUpTable, 
     *  sensorValues) 
	y(189)    = vGetSensorValue('YH21', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(190)    = vGetSensorValue('YH22', 
     *jSensorLookUpTable, 
     *sensorValues) 
	y(191)    = vGetSensorValue('YH23', 
     *jSensorLookUpTable, 
     *sensorValues) 
     	y(192)= vGetSensorValue('YH24', 
     *jSensorLookUpTable, 
     * sensorValues)   

           Cupx= vGetSensorValue('Z-10', 
     * jSensorLookUpTable,  
     * sensorValues)
          Corex= vGetSensorValue('Z-13', 
     * jSensorLookUpTable,  
     * sensorValues)
          Drivex= vGetSensorValue('Z-12', 
     * jSensorLookUpTable,  
     * sensorValues)
     !!     Cupy= vGetSensorValue('Cupy', 
     !!* jSensorLookUpTable,  
     !!* sensorValues)
     !!     Cdowny= vGetSensorValue('Cdowny', 
     !!* jSensorLookUpTable,  
     !!* sensorValues)

          
          
		R=1000.0*(Radius(x(1),x(9),x(17),y(1),y(9),y(17))
     * +Radius(x(2),x(10),x(18),y(2),y(10),y(18))
     * +Radius(x(3),x(11),x(19),y(3),y(11),y(19))
     * +Radius(x(4),x(12),x(20),y(4),y(12),y(20)) 
     * +Radius(x(5),x(13),x(21),y(5),y(13),y(21)) 
     * +Radius(x(6),x(14),x(22),y(6),y(14),y(22))  
     * +Radius(x(7),x(15),x(23),y(7),y(15),y(23)) 
     * +Radius(x(8),x(16),x(24),y(8),y(16),y(24))  
     * )/8
c 保存R1
          Write(700,*),time(iTotalTime),R
          open(700,File= 
     *'D:\abaquss\jilu2\
     *R.txt',status='unknown')

           R2=1000.0*(Radius(x(25),x(33),x(41),y(25),y(33),y(41))
     * +Radius(x(26),x(34),x(42),y(26),y(34),y(42))
     * +Radius(x(27),x(35),x(43),y(27),y(35),y(43))
     * +Radius(x(28),x(36),x(44),y(28),y(36),y(44)) 
     * +Radius(x(29),x(37),x(45),y(29),y(37),y(45)) 
     * +Radius(x(30),x(38),x(46),y(30),y(38),y(46))  
     * +Radius(x(31),x(39),x(47),y(31),y(39),y(47)) 
     * +Radius(x(32),x(40),x(48),y(32),y(40),y(48))  
     * )/8
c 保存R2
          Write(710,*),time(iTotalTime),R2
          open(710,File= 
     *'D:\abaquss\jilu2\
     *R2.txt',status='unknown')
          
          R3=1000.0*(Radius(x(49),x(57),x(65),y(49),y(57),y(65))
     * +Radius(x(50),x(58),x(66),y(50),y(58),y(66))
     * +Radius(x(51),x(59),x(67),y(51),y(59),y(67))
     * +Radius(x(52),x(60),x(68),y(52),y(60),y(68)) 
     * +Radius(x(53),x(61),x(69),y(53),y(61),y(69)) 
     * +Radius(x(54),x(62),x(70),y(54),y(62),y(70))  
     * +Radius(x(55),x(63),x(71),y(55),y(63),y(71)) 
     * +Radius(x(56),x(64),x(72),y(56),y(64),y(72))  
     * )/8
c 保存R3
          Write(720,*),time(iTotalTime),R3
          open(720,File= 
     *'D:\abaquss\jilu2\
     *R3.txt',status='unknown')                  

           R4=1000.0*(Radius(x(73),x(81),x(89),y(73),y(81),y(89))
     * +Radius(x(74),x(82),x(90),y(74),y(82),y(90))
     * +Radius(x(75),x(83),x(91),y(75),y(83),y(91))
     * +Radius(x(76),x(84),x(92),y(76),y(84),y(92)) 
     * +Radius(x(77),x(85),x(93),y(77),y(85),y(93)) 
     * +Radius(x(78),x(86),x(94),y(78),y(86),y(94))  
     * +Radius(x(79),x(87),x(95),y(79),y(87),y(95)) 
     * +Radius(x(80),x(88),x(96),y(80),y(88),y(96))  
     * )/8
c 保存R4
          Write(730,*),time(iTotalTime),R4
          open(730,File= 
     *'D:\abaquss\jilu2\
     *R4.txt',status='unknown')   
          
           R5=1000.0*(Radius(x(97),x(105),x(113),y(97),y(105),y(113))
     * +Radius(x(98),x(106),x(114),y(98),y(106),y(114))
     * +Radius(x(99),x(107),x(115),y(99),y(107),y(115))
     * +Radius(x(100),x(108),x(116),y(100),y(108),y(116)) 
     * +Radius(x(101),x(109),x(117),y(101),y(109),y(117)) 
     * +Radius(x(102),x(110),x(118),y(102),y(110),y(118))  
     * +Radius(x(103),x(111),x(119),y(103),y(111),y(119)) 
     * +Radius(x(104),x(112),x(120),y(104),y(112),y(120))  
     * )/8
c 保存R5
          Write(740,*),time(iTotalTime),R5
          open(740,File= 
     *'D:\abaquss\jilu2\
     *R5.txt',status='unknown')
          
          R6=1000.0*(Radius(x(121),x(129),x(137),y(121),y(129),y(137))
     * +Radius(x(122),x(130),x(138),y(122),y(130),y(138))
     * +Radius(x(123),x(131),x(139),y(123),y(131),y(139))
     * +Radius(x(124),x(132),x(140),y(124),y(132),y(140)) 
     * +Radius(x(125),x(133),x(141),y(125),y(133),y(141)) 
     * +Radius(x(126),x(134),x(142),y(126),y(134),y(142))  
     * +Radius(x(127),x(135),x(143),y(127),y(135),y(143)) 
     * +Radius(x(128),x(136),x(144),y(128),y(136),y(144))  
     * )/8
c 保存R6
          Write(750,*),time(iTotalTime),R6
          open(750,File= 
     *'D:\abaquss\jilu2\
     *R6.txt',status='unknown')
          
          R7=1000.0*(Radius(x(145),x(153),x(161),y(145),y(153),y(161))
     * +Radius(x(146),x(154),x(162),y(146),y(154),y(162))
     * +Radius(x(147),x(155),x(163),y(147),y(155),y(163))
     * +Radius(x(148),x(156),x(164),y(148),y(156),y(164)) 
     * +Radius(x(149),x(157),x(165),y(149),y(157),y(165)) 
     * +Radius(x(150),x(158),x(166),y(150),y(158),y(166))  
     * +Radius(x(151),x(159),x(167),y(151),y(159),y(167)) 
     * +Radius(x(152),x(160),x(168),y(152),y(160),y(168))  
     * )/8
c 保存R7
          Write(760,*),time(iTotalTime),R7
          open(760,File= 
     *'D:\abaquss\jilu2\
     *R7.txt',status='unknown')
          
          R8=1000.0*(Radius(x(169),x(177),x(185),y(169),y(177),y(185))
     * +Radius(x(170),x(178),x(186),y(170),y(178),y(186))
     * +Radius(x(171),x(179),x(187),y(171),y(179),y(187))
     * +Radius(x(172),x(180),x(188),y(172),y(180),y(188)) 
     * +Radius(x(173),x(181),x(189),y(173),y(181),y(189)) 
     * +Radius(x(174),x(182),x(190),y(174),y(182),y(190))  
     * +Radius(x(175),x(183),x(191),y(175),y(183),y(191)) 
     * +Radius(x(176),x(184),x(192),y(176),y(184),y(192))  
     * )/8
c 保存R8
          Write(770,*),time(iTotalTime),R8
          open(770,File= 
     *'D:\abaquss\jilu2\
     *R8.txt',status='unknown')
    
    !R2~R8 对应环件不同高度的 7 个截面，覆盖整个环件高度，用于监测环件外表面的不平整度      
          
c  计算测量凹槽半径的方差
          
      !环件椭圆度监测
          
     !! ssa=sqrt(((1000*Radius(x(1),x(2),x(3),y(1),y(2),y(3))-R)**2
     !!* +(1000*Radius(x(4),x(5),x(6),y(4),y(5),y(6))-R)**2
     !!* +(1000*Radius(x(7),x(8),x(9),y(7),y(8),y(9))-R)**2
     !!* +(1000*Radius(x(10),x(11),x(12),y(10),y(11),y(12))-R)**2 
     !!* +(1000*Radius(x(13),x(14),x(15),y(13),y(14),y(15))-R)**2 
     !!* +(1000*Radius(x(16),x(17),x(18),y(16),y(17),y(18))-R)**2  
     !!* +(1000*Radius(x(19),x(20),x(21),y(19),y(20),y(21))-R)**2 
     !!* +(1000*Radius(x(22),x(23),x(24),y(22),y(23),y(24))-R)**2)/8
     !!*)
     !!    Write(780,*),time(iTotalTime),ssa
     !!     open(780,File= 
     !!*'D:\abaquss\jilu2\
     !!*ssa.txt',status='unknown')
     !
!c  计算测量上凸台半径的方差       
!      sst=sqrt(((1000*Radius(x(31),x(32),x(33),y(31),y(32),y(33))-Rt)**2
!     * +(1000*Radius(x(34),x(35),x(36),y(34),y(35),y(36))-Rt)**2
!     * +(1000*Radius(x(37),x(38),x(39),y(37),y(38),y(39))-Rt)**2
!     * +(1000*Radius(x(40),x(41),x(42),y(40),y(41),y(42))-Rt)**2 
!     * +(1000*Radius(x(43),x(44),x(45),y(43),y(44),y(45))-Rt)**2 
!     * +(1000*Radius(x(46),x(47),x(48),y(46),y(47),y(48))-Rt)**2  
!     * +(1000*Radius(x(49),x(50),x(51),y(49),y(50),y(51))-Rt)**2 
!     * +(1000*Radius(x(52),x(53),x(54),y(52),y(53),y(54))-Rt)**2)/8
!     *)
!      
c  计算理想半径！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
      Rideal=(2*R0-b0)*b0*h0/2/(b*1000)/(h0*1000)
     *+500*b

      
           
c 导辊旋转中心到环件中心距离！！！！！！！！！！！！！！！！！！！！！！！！！          !用什么算的半径这里就用什么
		GRotRigD=(Rmain+R-GRotXo)**2+GRotYo**2
                    

c 导辊旋转角度！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！          !用什么算的半径这里就用什么
		w=atan((-R-Rmain+GRotXo)/GRotYo)
     *  +acos(0.5*(GRotRigD+Rg**2-(Rguide+R)**2)/(sqrt(GRotRigD)*Rg))

          
	angle=w 
        Write(1100,*),time(iTotalTime),angle,w
        open(1100,File=
     *'D:\abaquss\jilu2\
     *angle.txt',status='unknown')   

c 厚度测量值！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
      bdet=Corex-Drivex-Rmain/1000-Ridle/1000
      If(bdet>=bf/1000) then
          b=bdet
      Else
          b=bf/1000
      End If
      
          Write(1400,*),time(iTotalTime),b
          open(1400,File= 
     *'D:\abaquss\jilu2\
     *b.txt',status='unknown')
c 环件直径长大速度！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
      h=Cupy-Cdowny
      If(b>=0.261) then
                Vd0=0.8*time(iTotalTime)
          If( Vd0<=8.) Then
              Vd=Vd0
          Else 
          Vd=8.
              End If
      Else if (0.261>b.and.b>=0.181) then
      Vd=8.
      Else 
          Vd=0.4*(b*1000-161)
      End if
     
!c 轧制曲线！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
!                 
!       If(b>0.21875) then
!c 上凸型！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
!           n=2.5
!       Vf=Vd/((2*R0-b0)*b0*h0*(1/(b*1000.)**2/(h*1000.)+n*((h0-hi)/
!     *(b0-bi)**n)*(b0-b*1000)**(n-1)/(b*1000.)/(h*1000.)**2)-1)        !芯辊径向速度(mm)
!      Va=n*((h0-hi)/(b0-bi)**n)*(b0-b*1000)**(n-1)*Vf !轧制曲线求导
!c 下凹型！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
!      Else If (0.21875>=b.and.b>0.161) Then
!            n=2.5
!      Vf=Vd/((2*R0-b0)*b0*h0*(1/(b*1000.)**2/(h*1000.)+n*((hi-hf)/
!     *(bi-bf)**n)*(b*1000-bf)**(n-1)/(b*1000.)/(h*1000.)**2)-1)  
!      Va=n*((hi-hf)/(bi-bf)**n)*(b*1000-bf)**(n-1)*Vf 
!      Write(700,*) Vf
!       Else
!      Vf=0
!      Va=0
!       End if
c       Write(800,*),time(iTotalTime),b  
c       curve=h/b
c       Write(900,*),b,curve
c       Write(1100,*),time(iTotalTime),Va   
c	Write(1200,*),time(iTotalTime),Vf 
c 芯辊！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
						if(lFirst.eq.1) then
c First call of this channel: remember the current value in its own
c previous-value slot and output zero, so no artificial jump is applied.
							if(ampName(1:2).eq.'A1') then
								svars(5)=GRotYo+Rg*cos(angle)
							else if(ampName(1:2).eq.'A2') then
								svars(6)=GRotXo+Rg*sin(angle)
							else if(ampName(1:3).eq.'A3-') then
								svars(10)=2*R
							else if(ampName(1:2).eq.'A3') then
								svars(7)=2*R
							else if(ampName(1:3).eq.'A4-') then
								svars(11)=2*R
							else if(ampName(1:2).eq.'A4') then
								svars(8)=2*R
							end if
							ampValueNew=0.0
						else
							if(t.ne.0) then
							select case(ampName) 
c      case('CORE')
c          ampValueNew=Vf
c 导辊Z方向！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
	case('A1')
	svars(9)=GRotYo+Rg*cos(angle)       
	ampValueNew=(svars(9)-svars(5))/t 
	if(ampValueNew.gt.0.0) then 
	ampValueNew=0 
      end if
      if(time(iTotalTime).le.0.1) then
      ampValueNew=0
      end if

       svars(5)=svars(9)  
c 导辊X方向！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
	case('A2')                           
	svars(9)=GRotXo+Rg*sin(angle)       
	ampValueNew=(svars(9)-svars(6))/t 
      if(time(iTotalTime).le.0.1) then
      ampValueNew=0
      end if
	if(w.ge.pi/2) then

		if(ampValueNew.ge.0.0) then 
		ampValueNew=0 
		end if 
	else
	    if(ampValueNew.le.0.0) then 
		ampValueNew=0 
		end if 
      end if

	svars(6)=svars(9) 
      
      Write(1500,*),time(iTotalTime),svars(6)
          open(1500,File= 
     *'D:\abaquss\jilu2\
     *daogun.txt',status='unknown')        

          
          
          
c 上锥辊X！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！局部坐标系
	       case('A3') 
	svars(9)=2*R
	ampValueNew=(svars(9)-svars(7))*cos(0.279)
     *         /t!-Va *sin(Acone)
c 幅值保护：合成速度限15mm/s，A3=cos分量限14.42，后退方向同样保护（可调）
        if(ampValueNew.gt. 14.42) ampValueNew= 14.42
        if(ampValueNew.lt.-14.42) ampValueNew=-14.42
       Write(1600,*),'上一步直径=',svars(7),'本步直径=',
     *svars(9),'锥辊（上）X方向速度=（本步-上步）/t=',ampValueNew
								svars(7)=svars(9) 
                                  
                           
c 上锥辊Y！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！                                 
								case('A4') 
								svars(9)=2*R
								ampValueNew=(svars(9)-svars(8))*sin(0.279)
     *		/t!+Va*cos(Acone) 
c 幅值保护：A4=sin分量限4.13，与A3按16°比例同步限制，方向不变（可调）
        if(ampValueNew.gt. 4.13) ampValueNew= 4.13
        if(ampValueNew.lt.-4.13) ampValueNew=-4.13
                   Write(1700,*),'上一步直径=',svars(8),'本步直径=',
     *svars(9),'锥辊（上）Y方向速度=（本步-上步）/t=',ampValueNew
                   Write(1300,*),time(iTotalTime),(svars(9)-svars(8))/t
								svars(8)=svars(9)
          Write(1400,*),time(iTotalTime),ampValueNew
          open(1400,File= 
     *'D:\abaquss\jilu2\
     *sst1.txt',status='unknown') 
c 下锥辊X！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
      case('A3-1') 
	svars(9)=2*R
	ampValueNew=(svars(9)-svars(10))*cos(0.279)
     *         /t
c 幅值保护：合成速度限15mm/s，A3-1=cos分量限14.42，后退方向同样保护（可调）
        if(ampValueNew.gt. 14.42) ampValueNew= 14.42
        if(ampValueNew.lt.-14.42) ampValueNew=-14.42
      Write(1800,*),'上一步直径=',svars(10),'本步直径=',
     *svars(9),'锥辊（下）X方向速度=（本步-上步）/t=',ampValueNew
								svars(10)=svars(9) 
c 下锥辊Y！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！                                 
      case('A4-1') 
								svars(9)=2*R
								ampValueNew=(svars(9)-svars(11))*sin(0.279)
     *		/t 
c 幅值保护：A4-1=sin分量限4.13，与A3-1按16°比例同步限制，方向不变（可调）
        if(ampValueNew.gt. 4.13) ampValueNew= 4.13
        if(ampValueNew.lt.-4.13) ampValueNew=-4.13
          Write(1900,*), '上一步直径=',svars(11),'本步直径=',
     *svars(9),'锥辊（下）Y速度=（本步-上步）/t=',ampValueNew
								svars(11)=svars(9)
c 锥辊W！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！                                 
                                  case('A5')  
                          Lcone=2*R-1000.*Cupx+1000.*Corex-Ridle-b*1000.
                          L=(Lcone-b*730.)*sin(16*pi/180) 
                          ampValueNew=0.4*Vdrive*(R-b*730.)/R/L
							end select
							end if
						end if
				else 
					ampValueNew=ampValueOld 
      end if   
      end if
	return 
	end 


	FUNCTION Radius(x1,x2,x3,y1,y2,y3)                                                               
	real(8) Radius,x1,x2,x3,y1,y2,y3,xo,yo,t 
	t=x1*y2-x1*y3-y1*x2+y1*x3-y2*x3+y3*x2 
	if(t .ne. 0.0) then 
	xo=0.5*(-x2**2*y1+y3*x2**2+y2*y1**2-y2*x3**2+x3**2*y1+y3*y2**2- 
     *  y2**2*y1+y2*x1**2-y3*y1**2-y3*x1**2+y3**2*y1-y2*y3**2)/t             
	yo=-0.5*(-x1*x2**2-y2**2*x1+x1*x3**2+y3**2*x1+x1**2*x2-x1**2*x3+ 
     *      y1**2*x2-y1**2*x3+x2**2*x3+y2**2*x3-x2*x3**2-x2*y3**2)/t    


	Radius=sqrt((xo-x1)**2+(yo-y1)**2) 
	end if   
	END 
