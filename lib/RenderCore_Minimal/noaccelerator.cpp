bool NoAccelerator::Intersect( const Ray &r, HitInfo *hitinfo ) const
{
	throw runtime_error( "NoAccelerator::Intersect not implemented." );
	for ( const auto p : primitives )
	{
		p->Intersect( r, hitinfo );
	}
}

bool NoAccelerator::IntersectP( const Ray &r ) const
{
	throw runtime_error( "NoAccelerator::IntersectP not implemented." );
	for ( const auto p : primitives )
	{
		p->IntersectP(r);
	}
}