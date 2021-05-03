    // Output the region flag (0: recirculation; 1: contraction; 2: upper channel) 

    // Region Info for pehill case
    //forAll(mesh.C(), cellI)
    //{
    //    scalar x = mesh.C()[cellI].component(0);
    //    scalar y = mesh.C()[cellI].component(1);
    //    scalar z = mesh.C()[cellI].component(2);
    //    if (y <= 1.2 && x <= 5)
    //    {
    //        regionFlag[cellI] = 1;
    //    }
    //    else
    //    {
    //        regionFlag[cellI] = 0;
    //    }
    //}
    //regionFlag.write();
    
    // Region Info for channel case
    //forAll(mesh.C(), cellI)
    //{
    //    scalar y = mesh.C()[cellI].component(1);
    //    if (y < 5.0 / 5200.0)
    //    {
    //        regionFlag[cellI] = 0;
    //    }
    //    else if(y < 30.0 / 5200.0)
    //    {
    //        regionFlag[cellI] = 1;
    //    }
    //    else
    //    {
    //        regionFlag[cellI] = 2;
    //    }
    //}
    //regionFlag.write();

    // Region Info for wavy channel case
    //forAll(mesh.C(), cellI)
    //{
    //    scalar x = mesh.C()[cellI].component(0);
    //    scalar y = mesh.C()[cellI].component(1);
    //    scalar z = mesh.C()[cellI].component(2);
    //    if (y < 0.1 && x < 0.5 )
    //    {
    //        regionFlag[cellI] = 1;
    //    }
    //    else
    //    {
    //        regionFlag[cellI] = 0;
    //    }
    //}
    //regionFlag.write();

    // Region Info for backward-facing step case
    //forAll(mesh.C(), cellI)
    //{
    //    scalar x = mesh.C()[cellI].component(0);
    //    scalar y = mesh.C()[cellI].component(1);
    //    scalar z = mesh.C()[cellI].component(2);
    //    if (y > 1)
    //    {
    //        regionFlag[cellI] = 2;
    //    }
    //    else if(x < 10 )
    //    {
    //        regionFlag[cellI] = 0;
    //    }
    //    else
    //    {
    //        regionFlag[cellI] = 1;
    //    }
    //}
    //regionFlag.write();

    // Region Info for curved-backward-facing step case
    //forAll(mesh.C(), cellI)
    //{
    //    scalar x = mesh.C()[cellI].component(0);
    //    scalar y = mesh.C()[cellI].component(1);
    //    scalar z = mesh.C()[cellI].component(2);
    //    if (y <= 1.5 && x >= 0 && x <=5)
    //    {
    //        regionFlag[cellI] = 1;
    //    }
    //    else
    //    {
    //        regionFlag[cellI] = 0;
    //    }
    //}
    //regionFlag.write();

    // Region Info for curved-backward-facing step case
    forAll(mesh.C(), cellI)
    {
        scalar x = mesh.C()[cellI].component(0);
        scalar y = mesh.C()[cellI].component(1);
        scalar z = mesh.C()[cellI].component(2);
        if (y <= 1 && x >= 5 && x <=8)
        {
            regionFlag[cellI] = 1;
        }
        else
        {
            regionFlag[cellI] = 0;
        }
    }
    regionFlag.write();
