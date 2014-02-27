calcProfileWidth <- function(theta_a, theta_s, r){
        if(theta_a > 2*pi | theta_a < 0) 
		stop('theta_a is out of bounds. theta_a should be in 0<a<2*pi')
        if(theta_s > 2*pi | theta_s < 0) 
		stop('theta_s is out of bounds. theta_s should be in 0<a<2*pi')

	if(theta_a > pi){
	        if(theta_a < 4*pi - 2*theta_s){
		        p <- r*(theta_s - cos(theta_a/2) + 1)/pi
                } else if(theta_a <= 3*pi - theta_s){
                        p <- r*(theta_s - cos(theta_a/2) + cos(theta_a/2 + theta_s))/pi
                } else {
                        p <- r*(theta_s + 2*sin(theta_s/2))/pi
                }
        } else {
        	if(theta_a < 4*pi - 2*theta_s){
                        p <- r*(theta_s*sin(theta_a/2) - cos(theta_a/2) + 1)/pi
 		} else {
                        p <- r*(theta_s*sin(theta_a/2) - cos(theta_a/2) + cos(theta_a/2 + theta_s))/pi
                }
        }
        return(p)
}