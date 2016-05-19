int main() {
Accumulation fold (Accumulator f, Accumulation x0, Observations zs)
{   for (zs.current = 0; zs.current < zs.count; ++zs.current)
    {   x0 = f (x0, zs.observations[zs.current]);   }
    return x0;   }
return 0;
}
