function t = rostime2sec(stamp)
    t = stamp.Sec + 1e-9 * stamp.Nsec;
end