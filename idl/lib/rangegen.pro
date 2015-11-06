  function rangegen, range, n

; generates equidistant grid of n points from range(0) to range(1)

    return, findgen(n)/(n-1.)*(range(1)-range(0)) + range(0)

  end
