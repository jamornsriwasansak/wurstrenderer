#pragma once

template <typename Type>
struct ScopedAssignment
{
	// a very useful class that is copied from pbrt-v3
	ScopedAssignment(Type *target = nullptr, Type value = Type()) : target(target)
	{
		if (target)
		{
			backup = *target;
			*target = value;
		}
	}

	~ScopedAssignment()
	{
		if (target) *target = backup;
	}

	ScopedAssignment(const ScopedAssignment &) = delete;
	ScopedAssignment &operator=(const ScopedAssignment &) = delete;
	ScopedAssignment &operator=(ScopedAssignment &&other)
	{
		target = other.target;
		backup = other.backup;
		other.target = nullptr;
		return *this;
	}

	Type *target, backup;
};
