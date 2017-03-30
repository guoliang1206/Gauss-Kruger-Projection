

public class DmsRad {

	public static double p = 180.0 / Math.PI * 3600;
	
	// 将度分秒化为弧度值
	public static double dms2rad(int d, int m, int s){
		double dms = (d * 3600.0 + m * 60.0 + s) / p;
		return dms;
	}
	
	//将弧度值转化为度分秒
	public static int [] rad2dms (double rad){
		int []a = new int[3];
		double dms = rad *p;
		a[0] = (int) Math.floor(dms/3600.0);
		a[1] = (int) Math.floor((dms-a[0]*3600)/60.0);
		a[2] = (int) (dms-a[0]*3600)-a[1]*60;
		return a;
	}
	
	//将度分秒转化为十进制坐标
	public static double dms2dec (int [] dms){
		double dec = 0.0;
		dec = dms[0]+dms[1]/60.0+dms[2]/3600.0;
		return dec;
	}
	
	//将十进制坐标转换为度分秒
	public static int [] dec2Dms (double x){
		int []a = new int[3];
		
		double x1 = Math.abs(x);
		double degree,min,sec;
		degree = Math.floor(x1);
		x1 = x1-degree;
		min =  Math.floor(x1*60);
		x1 = (x1*60-min);
		sec = x1*60;
		a[0] = (int) degree;
		a[1] = (int) min;
		a[2] = (int) sec;
		
		return a;
	}
}
