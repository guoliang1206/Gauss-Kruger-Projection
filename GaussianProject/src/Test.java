public class Test {

	//中央子午线经度
	public static double L0 = 0;
	//WGS-84   椭球体参数
	public static double a = 6378137.0;	//major semi axis
	public static double efang = 0.0066943799901413;		//square of e
	public static double e2fang = 0.0067394967422764;	//suqre of e2
	

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		System.out.println("Test ");
		gaussPositiveFormula(36.082846593,119.3520122887);
	}

	/* 高斯投影正算
	 * params:B 纬度 ，L 经度 
	 * 十进制经纬度 转换为度分秒 度分秒转换为弧度
	 */
	public static void gaussPositiveFormula (double B,double L){
		
		//计算该地区的中央子午线经度
		L0 = Math.floor((L+1.5)/3)*3;
		double l =  L - L0;	//  P点所在经度L与中央子午线L0经度之差	

		//坐标位置转化为弧度
		int [] array_B = DmsRad.dec2Dms(B);
		B = DmsRad.dms2rad(array_B[0], array_B[1], array_B[2]);
		int [] array_L = DmsRad.dec2Dms(L);
		L = DmsRad.dms2rad(array_L[0], array_L[1], array_L[2]);
		int [] array_L0 = DmsRad.dec2Dms(L0);
		L0 = DmsRad.dms2rad(array_L0[0], array_L0[1], array_L0[2]);	
		
		System.out.println("B = "+ B);
		System.out.println("L = "+ L);
		
		
		//将l转换为以秒为单位
		l = formToDegress(l);
			
		double c = 6399593.6258;
		double sinB = Math.sin(B);
		double cosB = Math.cos(B);
		double t = Math.tan(B);
		double p = 180*3600/Math.PI;
		double nfang = e2fang*cosB*cosB;

		double V = Math.pow((1+e2fang*cosB*cosB), 1/2.0);
		double N = a*Math.pow((1-efang*sinB*sinB),-1/2);   
		N = c/V;

		//主曲率计算
		double m0,m2,m4,m6,m8;
		m0 = a*(1-efang);
		m2 = 3.0/2.0*efang*m0;
		m4 = efang*m2*5.0/4.0;
		m6 = efang*m4*7.0/6.0;
		m8 = efang*m6*9.0/8.0;

		//子午线曲率计算
		double a0,a2,a4,a6,a8;
		a0 = m0 + m2/2.0 + m4*3.0/8.0 + m6*5.0/16.0 + m8*35.0/128.0; 
		a2 = m2/2.0 + m4/2.0 + m6*15.0/32.0 + m8*7.0/16.0;
		a4 = m4/8.0 + m6*3.0/16.0 + m8*7.0/32.0;
		a6 = m6/32.0 + m8/16.0;
		a8 = m8/128.0;
	

		//中央子午线半径计算
		double X0 = a0*B -a2*Math.sin(2*B)/2.0 + a4*Math.sin(4*B)/4.0 - a6*Math.sin(6*B)/6.0 + a8*Math.sin(8*B)/8.0;
		System.out.println("X0 = "+ X0);

		double x = X0 + N*sinB*cosB*l*l/(p*p*2.0) + N*sinB*Math.pow(cosB, 3)*(5-t*t+9*nfang+4*nfang*nfang)*Math.pow(l, 4)/(24*Math.pow(p, 4)) + N*sinB*Math.pow(cosB, 5)*(61-58*t*t + Math.pow(t, 4))*Math.pow(l, 6)/(720*Math.pow(p, 6)); 
		System.out.println("x = "+ x);
		double y = N*cosB*l/p  +  N*Math.pow(cosB, 3)*(1-t*t + nfang)*Math.pow(l, 3)/(6*Math.pow(p, 3)) + N*Math.pow(cosB, 5)*(5-18*t*t +Math.pow(t, 4) + 14*nfang -58*nfang*t*t)*Math.pow(l, 5)/(120*Math.pow(p, 5));
		
		//避免负坐标，添加500Km偏移量
		y = y+500000;
		System.out.println("y = "+ y);
		System.out.println(N*cosB*l/p );
		
		gaussInverseFormula(x,y,L0);
		
	}
	
	//x = 3994922.0539079206 	y = 441635.7722969403
	
	/* 高斯投影反算
	 * params: x，y ，高斯平面坐标点
	 **/
	public static void gaussInverseFormula(double x,double y,double L0){
		
		y = y - 500000;
	
		//主曲率计算
		double m0,m2,m4,m6,m8;
		m0 = a*(1-efang);
		m2 = 3.0/2.0*efang*m0;
		m4 = efang*m2*5.0/4.0;
		m6 = efang*m4*7.0/6.0;
		m8 = efang*m6*9.0/8.0;

		//子午线曲率计算
		double a0,a2,a4,a6,a8;
		a0 = m0 + m2/2.0 + m4*3.0/8.0 + m6*5.0/16.0 + m8*35.0/128.0; 
		a2 = m2/2.0 + m4/2.0 + m6*15.0/32.0 + m8*7.0/16.0;
		a4 = m4/8.0 + m6*3.0/16.0 + m8*7.0/32.0;
		a6 = m6/32.0 + m8/16.0;
		a8 = m8/128.0;
		
		double X = x;
		double FBf = 0;
		double Bf0 = X/a0,Bf1 = 0;
		System.out.println("Initial Bf0 = " + Bf0);
		System.out.println("x = "+ x + "y = "+ y);
		
		
		//计算Bf的值，直到满足条件
		while ((Bf0-Bf1)>=0.0001){
			Bf1 = Bf0;
			FBf = -a2*Math.sin(2*Bf0)/2 +a4*Math.sin(4*Bf0)/4 -a6*Math.sin(6*Bf0)/6 + a8*Math.sin(8*Bf0)/8;
			Bf0 = (X-FBf)/a0;
		}
		
		double Bf  = Bf0;
		//计算公式中参数
		double Wf = Math.sqrt(1-efang*Math.sin(Bf)*Math.sin(Bf));
		double Nf = a/Wf;
		double Mf = a*(1-efang)/Math.pow(Wf, 3);
		double nffang = e2fang*Math.cos(Bf)*Math.cos(Bf);
		double tf = Math.tan(Bf);
	
		double B = Bf - tf*y*y/(2*Mf*Nf) + tf*(5+3*tf*tf+nffang - 9*nffang*tf*tf)*Math.pow(y, 4)/(24*Mf*Math.pow(Nf, 3))- tf*(61+90*tf*tf+45*Math.pow(tf, 4))*Math.pow(y, 6)/(720*Mf*Math.pow(Nf, 5));
		
		double l = y/(Nf*Math.cos(Bf)) - (1+2*tf*tf+nffang)*Math.pow(y, 3)/(6*Math.pow(Nf, 3)*Math.cos(Bf)) + (5+28*tf*tf +24*Math.pow(tf, 4))*Math.pow(y, 5)/(120*Math.pow(Nf, 5)*Math.cos(Bf));
		
		double L =l+L0;
		
		System.out.println("L0 = "+ L0);
		System.out.println("B = 0.6297632754876651 =                     "+ B);
		System.out.println("L = 2.0830843992129098 =                     "+ L);
		
		//转化成为十进制经纬度格式
		int [] array_B = DmsRad.rad2dms(B);
		int [] array_L =DmsRad.rad2dms(L);
		double Bdec = DmsRad.dms2dec(array_B);
		double Ldec = DmsRad.dms2dec(array_L);
		
		System.out.println("B = "+ Bdec + "\nL = "+ Ldec);
	}


	/*
	 * 将坐标点转化为度数（秒）
	 * */
	public static double formToDegress (double x){

		double x1 = Math.abs(x);
		x1 = x;
		double degree,min,sec;

		degree = Math.floor(x1);
		x1 = x1-degree;
		min =  Math.floor(x1*60);
		x1 = (x1*60-min);
		sec = x1*60;

		double s = degree*3600+min*60+sec;

		return s;
	}

}
