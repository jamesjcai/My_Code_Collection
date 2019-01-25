function [D]= random_walk_3d_simulation ( step_num, walk_num )

%*****************************************************************************80
%
%% RANDOM_WALK_3D_SIMULATION simulates a random walk in 3D.
%
%  Discussion:
%
%    The expectation should be that, the average distance squared D^2 
%    is equal to the time, or number of steps N.
%
%    Or, equivalently
%
%      average ( D ) = sqrt ( N )
%
%    The program makes a plot of both the average and the maximum values
%    of D^2 versus time.  The maximum value grows much more quickly,
%    and that curve is much more jagged.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    07 February 2014
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer STEP_NUM, the number of steps to take in one test.
%
%    Input, integer WALK_NUM, the number of times to repeat the walk.
%
  if ( nargin < 1 )
    step_num = 500;
    fprintf ( 1, '\n' );
    fprintf ( 1, '  STEP_NUM not supplied, using default value %d\n', step_num );
  end

  if ( nargin < 2 )
    walk_num = 10;
    fprintf ( 1, '\n' );
    fprintf ( 1, '  WALK_NUM not supplied, using default value %d\n', walk_num );
  end
%
%  Set up arrays for plotting.
%
  time = 0 : step_num;
  d2_ave = zeros(step_num+1,1);
  d2_max = zeros(step_num+1,1);
%
%  Take the walk WALK_NUM times.
%
  for walk = 1 : walk_num

    x = zeros(step_num+1,1);
    y = zeros(step_num+1,1);
    z = zeros(step_num+1,1);
  
   for step = 2 : step_num + 1
%
%  We are currently at ( X(STEP-1), Y(STEP-1), Z(STEP-1) ).
%  Consider the six possible destinations.
%
      destination = [ x(step-1) + 1.0, y(step-1),       z(step-1); ...
                      x(step-1) - 1.0, y(step-1),       z(step-1); ...
                      x(step-1),       y(step-1) + 1.0, z(step-1); ...
                      x(step-1),       y(step-1) - 1.0, z(step-1); ...
                      x(step-1),       y(step-1),       z(step-1) + 1.0; ...
                      x(step-1),       y(step-1),       z(step-1) - 1.0 ];
%
%  Choose a destination.
%
      k = ceil ( 6.0 * rand ( 1, 1 ) );
%
%  Move there.
%
      x(step) = destination(k,1);
      y(step) = destination(k,2);
      z(step) = destination(k,3);
%
%  Update the sum of every particle's distance at step J.
%
      d2 = x(step)^2 + y(step)^2 + z(step)^2;
      d2_ave(step) = d2_ave(step) + d2;
      d2_max(step) = max ( d2_max(step), d2 );
 
   end

  end
  C=[x,y,z];D=C';
end