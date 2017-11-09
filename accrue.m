function Y = accrue(X,dir)
  % ACCRUE: function to accrue the values in a vector
  %   ACCRUE(X,dri):
  %   The first index will always be 0 and the second the first values
  %   Outputs the sum of all previous values in a vector.

  m = length(X)+1;
  if dir == 1
    Y=zeros(m,1);
    for n = 1:length(X)
      Y(n+1) = Y(n)+X(n);
    end
  end

  if dir == -1

    Y=zeros(m,1);
    for n = length(X):-1:1
      Y(m-n+1) = Y(m-n)+X(n);
    end
  end
