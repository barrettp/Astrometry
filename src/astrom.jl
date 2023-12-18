#   Astrometry

function aberration(target, observer, distance, lorenz)
    pv = target.position*observer.velocity
    w1, w2 = 1.0 + pv/(1.0 + lorenz), SCHWARZRADIUS/distance
    pp = [pos*lorenz + w1*vel + w2*(vel - pv*pos)
         for (pos, vel) in zip(target.position, observer.velocity)]
    rr = sqrt(sum(pp.*pp))
    pp./rr
end
