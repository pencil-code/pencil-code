// Get divergence of velocity.
divu_shock()
{
    // Discard values which do not contain negative divergence.
    divu = divergence(UU)
    if divu < 0.0 {
        return -divu
    }
    else {
        return 0.0
    }
}

// Calculate local maximum.
max5_shock()
{
    return max5(SHOCK)
}
smooth_shock()
{
}

