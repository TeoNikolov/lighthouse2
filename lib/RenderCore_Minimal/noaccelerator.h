#pragma once

namespace lh2core
{

class NoAccelerator : public Aggregate
{
  public:
	NoAccelerator( std::vector<std::shared_ptr<Primitive>> &primitives ) : primitives( primitives ) {}
	virtual bool Intersect( const Ray &r, HitInfo *hitinfo ) const;
	virtual bool IntersectP( const Ray &r ) const;

  private:
	std::vector<std::shared_ptr<Primitive>> primitives;
};

} // namespace lh2core