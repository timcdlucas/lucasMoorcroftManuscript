
# Functions to calculate density.
#
# Tim C.D. Lucas, Elizabeth Moorcroft, Robin Freeman, Marcus J. Rowcliffe, Kate E. Jones.
#
# calcDensity is the main function to calculate density.
# It takes parameters z, alpha, theta, r, animalSpeed, t
# z - The number of camera/acoustic counts or captures.
# alpha - Call width in radians.
# theta - Sensor width in radians.
# r - Sensor range in metres.
# animalSpeed - Average animal speed in metres per second.
# t - Length of survey in sensor seconds i.e. number of sensors x survey duration.
#
# calcAbundance calculates abundance rather than density and requires an extra parameter
# area - In metres squared. The size of the region being examined.

using Base.Test


# Internal function to calculate profile width as described in the text

function calcProfileWidth(alpha, theta, r)
        @test alpha <=2*pi
        @test alpha > 0
        @test theta <=2*pi
        @test theta > 0

	if alpha > pi
	        if alpha < 4*pi - 2*theta
		        p = r*(theta - cos(alpha/2) + 1)/pi
                elseif alpha <= 3*pi - theta
                        p = r*(theta - cos(alpha/2) + cos(alpha/2 + theta))/pi
                else
                        p = r*(theta + 2*sin(theta/2))/pi
                end
        else
        	if alpha < 4*pi - 2*theta
                        p = r*(theta*sin(alpha/2) - cos(alpha/2) + 1)/pi
 		else
                        p = r*(theta*sin(alpha/2) - cos(alpha/2) + cos(alpha/2 + theta))/pi
                end
        end
        return p 
end

# Calculate a population density. See above for units etc.

function calcDensity(z, alpha, theta, r, animalSpeed, t)
        @test isa(z, Real)
        @test z > 0
        @test animalSpeed > 0
        @test isa(animalSpeed, Real)
        @test t > 0
        @test isa(t, Real)

        p = calcProfileWidth(alpha, theta, r)
        # if(p <= 0) stop('Calculated profile width is 0. We would therefore expect 0 captures. If z is not zero, then the density is undefined.')
        D = z/{animalSpeed*t*p}
        return D
end


# Calculate abundance rather than density.

function calcAbundance(z, alpha, theta, r, animalSpeed, t, area)
         @test area >= 0
         @test isa(area, Real)
         D = calcDensity(z, alpha, theta, r, animalSpeed, t)
         A = D*area
         return A
                        
end

